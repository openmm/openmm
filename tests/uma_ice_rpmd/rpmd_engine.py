"""
rpmd_engine.py – Pure-Python ring-polymer MD with batched UMA force evaluation.

No MPI, no LAMMPS, no i-PI server.  A single Python process holds all bead
positions, calls the UMA predictor once per time step (batched over all beads
in one GPU forward pass), and integrates with velocity-Verlet + PILE-L
thermostat.

Unit conventions (used throughout this module)
----------------------------------------------
  Positions    : Å
  Velocities   : Å fs⁻¹
  Forces       : eV Å⁻¹
  Mass         : amu (g mol⁻¹)
  Time step    : fs
  Energy       : eV
  Temperature  : K

Newton's 2nd law in these units:
  a [Å fs⁻²] = F [eV Å⁻¹] / m [amu]  ×  FORCE_CONV

  FORCE_CONV = 9.6485e-3
  Derivation:
    a [m/s²] = F [eV/Å] × 1.60218e-9 N/(eV/Å) / (m [amu] × 1.66054e-27 kg/amu)
             = F/m × 9.6485e17 m/s²
    → ×1e-10 [Å/m] × 1e-30 [s²/fs²]  →  9.6485e-3 Å fs⁻² per (eV Å⁻¹ amu⁻¹)

Kinetic energy:
  KE [eV] = ½ Σ m [amu] v² [Å² fs⁻²] / FORCE_CONV

Maxwell-Boltzmann velocity σ per component per atom:
  σ [Å fs⁻¹] = sqrt(k_B T × FORCE_CONV / m)

PILE-L references
-----------------
  Ceriotti, Parrinello, Markland, Manolopoulos,
  J. Chem. Phys. 133, 124104 (2010).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

log = logging.getLogger(__name__)

# ── Physical constants ─────────────────────────────────────────────────────────
KB         = 8.617333262e-5   # eV K⁻¹
HBAR       = 0.6582119514     # eV fs  (ℏ)
FORCE_CONV = 9.6485334e-3     # [eV Å⁻¹ amu⁻¹] → [Å fs⁻²]


# ── Configuration ──────────────────────────────────────────────────────────────
@dataclass
class RPMDConfig:
    """Parameters for the RPMD simulation."""
    n_beads:          int            = 32
    temperature:      float          = 243.0   # K
    dt:               float          = 0.1     # fs
    charge:           int            = 0
    spin:             int            = 1
    task_name:        str            = "omol"
    model_name:       str            = "uma-s-1p1"
    device:           str            = "cuda"
    seed:             Optional[int]  = 42
    # Centroid Langevin friction γ₀ [fs⁻¹].  Default 1e-3 fs⁻¹ = 1 ps⁻¹.
    gamma_centroid:   float          = 1.0e-3
    # Split bead batch into chunks of this size to cap GPU memory.
    # None = no splitting (full batch in one call).
    chunk_size:       Optional[int]  = None


# ── Engine ─────────────────────────────────────────────────────────────────────
class RPMDEngine:
    """
    Velocity-Verlet + PILE-L ring-polymer MD.

    Ring-polymer Hamiltonian
    ~~~~~~~~~~~~~~~~~~~~~~~~
    H = Σᵢ [ pᵢ²/(2m) + V(qᵢ) ] + ½ Σᵢ m ω_N² (qᵢ − q_{i+1})²

    where  ω_N = N k_B T / ħ  (ring-polymer bending frequency, rad fs⁻¹).

    Integration scheme (BAB + O)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    One step:
      B) half-velocity kick  v += (F_phys + F_spring) / m × FORCE_CONV × dt/2
      A) full position drift  q += v × dt  (+ molecular PBC wrap)
      B) recompute forces (batched UMA + ring springs)
         half-velocity kick  v += (F_phys + F_spring) / m × FORCE_CONV × dt/2
      O) PILE-L Ornstein-Uhlenbeck in normal-mode velocity space

    PILE-L thermostat
    ~~~~~~~~~~~~~~~~~
    Normal-mode frequencies:  ω_k = 2 ω_N sin(kπ/N),  k = 0 … N−1
    Friction per mode:
        k = 0  (centroid):  γ₀ = gamma_centroid [user-set]
        k > 0  (internal):  γ_k = 2 ω_k  (critical damping)
    O-U update:  u_k → c₁_k u_k + c₂_{k,a} ξ  where ξ ~ N(0,1)
        c₁_k = exp(−γ_k dt)
        c₂_{k,a} = sqrt((1 − c₁_k²) × k_B T × FORCE_CONV / m_a)
    """

    def __init__(self, atoms, config: RPMDConfig) -> None:
        """
        Parameters
        ----------
        atoms  : ase.Atoms
            Initial configuration (single bead; replicated internally).
        config : RPMDConfig
        """
        self.cfg  = config
        self.rng  = np.random.default_rng(config.seed)

        n = config.n_beads
        self.n_beads = n
        self.symbols = list(atoms.get_chemical_symbols())
        self.n_atoms = len(self.symbols)
        self.cell    = np.asarray(atoms.get_cell(), dtype=float).copy()   # (3,3) Å
        self.pbc     = np.asarray(atoms.get_pbc()).copy()
        self.masses  = atoms.get_masses().copy()                          # (n_atoms,) amu

        # Broadcast shapes for vectorised operations
        self._masses3 = self.masses[None, :, None]    # (1, n_atoms, 1)

        # ── Ring-polymer state arrays ──────────────────────────────────────────
        pos0     = np.asarray(atoms.get_positions(), dtype=float)
        self.pos = np.tile(pos0, (n, 1, 1)).copy()        # (n_beads, n_atoms, 3) Å
        self.vel = self._maxwell_boltzmann_velocities()    # (n_beads, n_atoms, 3) Å/fs
        self.forces = np.zeros_like(self.pos)              # eV/Å  (filled before first step)

        # ── Ring-polymer frequency ─────────────────────────────────────────────
        self.omega_N: float = n * KB * config.temperature / HBAR   # rad fs⁻¹

        # ── Normal-mode transform and PILE-L coefficients ─────────────────────
        self._C           = self._build_nm_matrix()            # (n, n) orthogonal
        self._c1, self._c2 = self._pile_l_coefficients()       # (n,) and (n, n_atoms)

        # ── FairChem predictor (lazy-loaded on first force call) ──────────────
        self._predict_unit  = None
        self._adata_cache: list = []
        self._dataset_name: str | None = None

        self.step_count: int = 0
        self._pe_ev:     float = 0.0     # most recent physical potential energy [eV]

        log.info(
            "RPMDEngine: %d beads × %d atoms, T=%.1f K, "
            "ω_N=%.5f rad/fs, dt=%.3f fs",
            n, self.n_atoms, config.temperature, self.omega_N, config.dt,
        )

    # ── Public API ─────────────────────────────────────────────────────────────

    def initialise_forces(self) -> None:
        """
        Evaluate forces at the initial positions.
        Must be called once before the first :meth:`step`.
        """
        F_phys, E = self._batched_uma_forces()
        F_spring  = self._ring_spring_forces()
        self.forces  = F_phys + F_spring
        self._pe_ev  = E
        log.info(
            "Initial forces computed: E_phys=%.6f eV, T_inst=%.1f K",
            E, self.kinetic_temperature_k,
        )

    def step(self) -> None:
        """Advance the simulation by one time step (dt fs)."""
        dt   = self.cfg.dt
        inv_m = 1.0 / self._masses3    # (1, n_atoms, 1) amu⁻¹

        # B) half-kick
        self.vel += self.forces * inv_m * FORCE_CONV * (dt * 0.5)

        # A) full drift
        self.pos += self.vel * dt
        self._molecular_pbc_wrap()

        # B) recompute forces + half-kick
        F_phys, E   = self._batched_uma_forces()
        F_spring     = self._ring_spring_forces()
        self.forces  = F_phys + F_spring
        self._pe_ev  = E
        self.vel    += self.forces * inv_m * FORCE_CONV * (dt * 0.5)

        # O) PILE-L thermostat
        self._pile_l_step()

        self.step_count += 1

    # ── Thermodynamic observables ──────────────────────────────────────────────

    @property
    def kinetic_energy_ev(self) -> float:
        """Classical KE = ½ Σ_{beads,atoms} m v²  [eV]."""
        return float(
            0.5 * np.sum(self._masses3 * self.vel ** 2) / FORCE_CONV
        )

    @property
    def potential_energy_ev(self) -> float:
        """Sum of physical (ML) potential energies over all beads [eV]."""
        return self._pe_ev

    @property
    def ring_spring_energy_ev(self) -> float:
        """Ring-polymer spring potential energy [eV]."""
        diff = self.pos - np.roll(self.pos, 1, axis=0)   # (n, n_atoms, 3) Å
        if np.any(self.pbc):
            diff = self._minimum_image(diff)
        # V_spring = ½ m ω_N² Δq²  →  [eV] = amu × (rad/fs)² × Å² / FORCE_CONV
        return float(
            0.5 * np.sum(self._masses3 * diff ** 2) * self.omega_N ** 2 / FORCE_CONV
        )

    @property
    def kinetic_temperature_k(self) -> float:
        """Instantaneous kinetic temperature from classical KE [K]."""
        n_dof = 3 * self.n_beads * self.n_atoms - 3   # subtract COM translation
        return 2.0 * self.kinetic_energy_ev / (n_dof * KB)

    @property
    def centroid_positions(self) -> np.ndarray:
        """Bead-averaged (centroid) positions, shape (n_atoms, 3) Å."""
        return self.pos.mean(axis=0)

    # ── Force evaluation ───────────────────────────────────────────────────────

    def _batched_uma_forces(self) -> tuple[np.ndarray, float]:
        """
        Evaluate physical (ML) forces for all beads in one batched UMA call.

        Returns
        -------
        forces : ndarray, shape (n_beads, n_atoms, 3), eV Å⁻¹
        total_energy : float, eV (sum over all beads)
        """
        import torch
        from fairchem.core.datasets.atomic_data import (
            AtomicData,
            atomicdata_list_to_batch,
        )
        from ase import Atoms as _Atoms

        if self._predict_unit is None:
            self._load_model()

        cfg  = self.cfg
        n    = self.n_beads
        task = self._dataset_name
        pbc  = bool(np.any(self.pbc))
        cell = self.cell if pbc else None

        # Build AtomicData templates on first call; update positions on subsequent calls
        if len(self._adata_cache) < n:
            self._adata_cache = []
            for b in range(n):
                pos_w = self._wrap_bead_positions(
                    np.ascontiguousarray(self.pos[b], dtype=np.float64)
                )
                at = _Atoms(symbols=self.symbols, positions=pos_w, pbc=pbc)
                at.info["charge"] = cfg.charge
                at.info["spin"]   = cfg.spin
                if cell is not None:
                    at.set_cell(cell)
                data = AtomicData.from_ase(
                    at, task_name=task, r_edges=False,
                    r_data_keys=["spin", "charge"],
                )
                data.sid = [f"bead-{b}"]
                self._adata_cache.append(data)
        else:
            for b in range(n):
                pos_w = self._wrap_bead_positions(
                    np.ascontiguousarray(self.pos[b], dtype=np.float64)
                )
                self._adata_cache[b].pos = torch.from_numpy(
                    np.ascontiguousarray(pos_w, dtype=np.float32)
                )
                if cell is not None:
                    self._adata_cache[b].cell = torch.from_numpy(
                        np.ascontiguousarray(cell, dtype=np.float32)
                    ).unsqueeze(0)

        # Optional chunked batching to cap peak GPU memory.
        # Example: RPMDConfig(chunk_size=8) for 32 beads → 4 chunks of 8.
        chunk = self.cfg.chunk_size or n
        all_forces   = np.zeros((n, self.n_atoms, 3), dtype=np.float64)
        total_energy = 0.0

        for start in range(0, n, chunk):
            end      = min(start + chunk, n)
            sub_list = self._adata_cache[start:end]
            n_sub    = end - start

            batch = atomicdata_list_to_batch(sub_list)
            pred  = self._predict_unit.predict(batch)

            e_np = pred["energy"].detach().cpu().numpy()   # (n_sub,)
            f_np = pred["forces"].detach().cpu().numpy()   # (n_sub × n_atoms, 3)
            bi   = batch.batch.numpy()                     # atom index → bead index

            total_energy += float(np.sum(e_np))
            for local_b in range(n_sub):
                all_forces[start + local_b] = f_np[bi == local_b]

        return all_forces, total_energy

    def _ring_spring_forces(self) -> np.ndarray:
        """
        Cartesian ring-polymer spring forces (eV Å⁻¹).

        F_spring[i] = m ω_N² (q_{i−1} + q_{i+1} − 2 q_i) / FORCE_CONV

        Derivation:
          Spring acceleration: a[i] = ω_N² Δq  [Å fs⁻²]
          Newton:              F [eV Å⁻¹] = m [amu] × a [Å fs⁻²] / FORCE_CONV
          Combined:            F = m ω_N² Δq / FORCE_CONV
        """
        q_prev = np.roll(self.pos,  1, axis=0)      # q_{i-1}
        q_next = np.roll(self.pos, -1, axis=0)      # q_{i+1}
        delta  = q_prev + q_next - 2.0 * self.pos   # (n_beads, n_atoms, 3) Å

        if np.any(self.pbc):
            delta = self._minimum_image(delta)

        return self._masses3 * self.omega_N ** 2 * delta / FORCE_CONV

    # ── PILE-L thermostat ──────────────────────────────────────────────────────

    def _build_nm_matrix(self) -> np.ndarray:
        """
        Real orthogonal normal-mode (NM) transformation matrix C[k, j], (N, N).

        Transforms bead coordinates to normal-mode coordinates:
            u[k] = Σ_j C[k,j] x[j]

        Convention (Ceriotti et al. 2010 / i-PI):
            k = 0:          C[0,j]    = 1/√N  (centroid)
            k = 1…N/2−1:    C[k,j]    = √(2/N) cos(2πkj/N)
            k = N/2 (even): C[N/2,j]  = (−1)ʲ/√N
            k > N/2:        C[k,j]    = √(2/N) sin(2π(N−k)j/N)
        """
        N = self.n_beads
        C = np.zeros((N, N))
        j = np.arange(N)

        C[0] = 1.0 / np.sqrt(N)

        for k in range(1, N // 2):
            C[k] = np.sqrt(2.0 / N) * np.cos(2.0 * np.pi * k * j / N)

        if N % 2 == 0:
            C[N // 2] = ((-1.0) ** j) / np.sqrt(N)

        for k in range(N // 2 + 1, N):
            C[k] = np.sqrt(2.0 / N) * np.sin(2.0 * np.pi * (N - k) * j / N)

        return C

    def _pile_l_coefficients(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Precompute PILE-L O-U coefficients for each normal mode.

        Returns
        -------
        c1 : ndarray, (n_beads,)
            Deterministic factor: c1[k] = exp(−γ_k dt)
        c2 : ndarray, (n_beads, n_atoms)
            Stochastic amplitude per mode per atom:
            c2[k,a] = sqrt((1 − c1[k]²) × k_B T × FORCE_CONV / m_a)
        """
        N  = self.n_beads
        dt = self.cfg.dt
        T  = self.cfg.temperature
        k  = np.arange(N)

        omega_k   = 2.0 * self.omega_N * np.sin(k * np.pi / N)  # (N,)  k=0 → 0
        gamma     = 2.0 * omega_k                                # (N,)  PILE-L critical friction
        gamma[0]  = self.cfg.gamma_centroid                      # centroid: user-set

        c1 = np.exp(-gamma * dt)                                 # (N,)

        # Equilibrium velocity σ per component per atom: sqrt(k_B T FORCE_CONV / m)
        sigma = np.sqrt(KB * T * FORCE_CONV / self.masses)       # (n_atoms,)
        c2    = np.sqrt(1.0 - c1 ** 2)[:, None] * sigma[None, :] # (N, n_atoms)

        return c1, c2

    def _pile_l_step(self) -> None:
        """
        Apply PILE-L O-U thermostat to velocities in normal-mode space.

        Algorithm (per step):
          1. Transform bead velocities → NM:   u[k] = Σ_j C[k,j] v[j]
          2. O-U update per mode per atom:      u[k] → c1[k] u[k] + c2[k,a] ξ
          3. Back-transform NM → bead:          v[j] = Σ_k C[k,j] u[k]  (C orthogonal)
        """
        # (n_beads, n_atoms, 3)  →  NM space
        u = np.einsum("kj,jad->kad", self._C, self.vel)

        # Stochastic update: broadcast c1 over (atoms, xyz) and c2 over xyz
        xi = self.rng.standard_normal(u.shape)
        u  = self._c1[:, None, None] * u + self._c2[:, :, None] * xi

        # Back-transform (C orthogonal → C⁻¹ = Cᵀ)
        self.vel = np.einsum("kj,kad->jad", self._C, u)

    # ── PBC utilities ──────────────────────────────────────────────────────────

    def _wrap_bead_positions(self, pos_ang: np.ndarray) -> np.ndarray:
        """
        Molecular PBC wrap for one bead's positions, shape (n_atoms, 3) Å.

        For water (OHH repeating):  wrap O into the primary cell, then place
        each H relative to its own O using the minimum-image vector.
        This preserves O−H bonds when the box is smaller than a bond length.

        Falls back to per-atom wrap for non-water systems.
        """
        if not np.any(self.pbc):
            return pos_ang

        n    = self.n_atoms
        cell = self.cell
        inv  = np.linalg.inv(cell)
        out  = np.array(pos_ang, dtype=np.float64, copy=True)

        is_water = (
            n % 3 == 0
            and all(
                self.symbols[3 * m]     == "O"
                and self.symbols[3 * m + 1] == "H"
                and self.symbols[3 * m + 2] == "H"
                for m in range(n // 3)
            )
        )

        if is_water:
            for m in range(n // 3):
                o   = out[3 * m].copy()
                fo  = o @ inv.T
                fo -= np.floor(fo)
                o_n = fo @ cell
                out[3 * m] = o_n
                for k in (1, 2):
                    d  = out[3 * m + k] - o
                    df = d @ inv.T
                    df -= np.round(df)
                    out[3 * m + k] = o_n + df @ cell
        else:
            for i in range(n):
                frac    = out[i] @ inv.T
                frac   -= np.floor(frac)
                out[i]  = frac @ cell

        return out

    def _molecular_pbc_wrap(self) -> None:
        """Wrap all bead positions in place (molecular wrap for water)."""
        if not np.any(self.pbc):
            return
        for b in range(self.n_beads):
            self.pos[b] = self._wrap_bead_positions(self.pos[b])

    def _minimum_image(self, delta: np.ndarray) -> np.ndarray:
        """
        Apply minimum-image convention to ring-polymer difference vectors.

        Parameters
        ----------
        delta : ndarray, (n_beads, n_atoms, 3) Å
        """
        if not np.any(self.pbc):
            return delta
        inv  = np.linalg.inv(self.cell)
        frac = delta @ inv.T
        frac -= np.round(frac)
        return frac @ self.cell

    # ── Velocity initialization ────────────────────────────────────────────────

    def _maxwell_boltzmann_velocities(self) -> np.ndarray:
        """
        Draw velocities from Maxwell-Boltzmann at config.temperature.

        Each bead is sampled independently.  Net linear momentum of the
        centroid is removed to prevent drift.

        Returns
        -------
        vel : ndarray, (n_beads, n_atoms, 3) Å fs⁻¹
        """
        T   = self.cfg.temperature
        vel = np.zeros((self.n_beads, self.n_atoms, 3))
        for a_idx, mass in enumerate(self.masses):
            # σ per Cartesian component: sqrt(k_B T FORCE_CONV / m)
            sigma = np.sqrt(KB * T * FORCE_CONV / mass)
            vel[:, a_idx, :] = self.rng.standard_normal((self.n_beads, 3)) * sigma

        # Remove net centroid linear momentum
        vel -= vel.mean(axis=(0, 1), keepdims=True)
        return vel

    # ── Model loading ──────────────────────────────────────────────────────────

    def _load_model(self) -> None:
        """Lazy-load the UMA predictor on the first force evaluation."""
        import torch
        from fairchem.core.calculate import pretrained_mlip

        cfg = self.cfg
        log.info("Loading UMA model '%s' on %s …", cfg.model_name, cfg.device)

        if cfg.device == "cuda":
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True

        predict_unit = pretrained_mlip.get_predict_unit(
            cfg.model_name,
            device=cfg.device,
        )

        valid = list(predict_unit.dataset_to_tasks.keys())
        if cfg.task_name in valid:
            self._dataset_name = cfg.task_name
        elif len(valid) == 1:
            self._dataset_name = valid[0]
        else:
            raise ValueError(
                f"task_name='{cfg.task_name}' not found; valid: {valid}"
            )

        self._predict_unit = predict_unit
        log.info("UMA ready (task: %s)", self._dataset_name)
