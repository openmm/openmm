# Fragmentation Techniques for Feasible ORCA Calculations and Force-Field Fitting

**Context:** Full MD box (e.g. protein + solvent) is too large for direct ORCA QM and for force-field optimization. This note summarizes fragmentation strategies from the literature and how to transfer them to a pipeline: **OpenMM MD → active learning → ORCA (QM) → ForceBalance/lightweight optimizer → updated force field**.

---

## 1. Why fragment?

- **ORCA** (and most QM codes) scale poorly with system size; full protein + solvent is intractable at DFT level.
- **ForceBalance**-style fitting with 4M structures is already heavy; using *fragment* QM data keeps reference calculations feasible and can still inform system-wide or residue-type parameters.
- **Goal:** Run QM only on *fragments* (e.g. &lt;150 atoms), then either (i) reconstruct full-system observables from fragment results, or (ii) fit force-field parameters to fragment-level reference data so they transfer to the full system.

---

## 2. Key fragmentation paradigms (from literature)

### 2.1 Automated QM region construction (Brunken & Reiher, *J. Chem. Theory Comput.* 2021)

**Paper:** Brunken, C.; Reiher, M. “Automated Construction of Quantum-Classical Hybrid Models.”  
**arXiv:** 2102.09355 | **DOI:** 10.1021/acs.jctc.1c00178

- **Idea:** Build **QM regions as spheres** around atoms (radius ~5.5–7 Å). All atoms inside the sphere form a fragment; **cut bonds are valence-saturated** (e.g. with H). Fragment sizes stay **under ~150 atoms** so DFT is feasible.
- **Cleavable bonds:** For biomolecules, C(sp³)–C and C(sp³)–N are natural cut sites. Rule-based saturation can be replaced later by more advanced embedding (e.g. frozen density).
- **Stochastic boundary:** When a cleavable bond is reached, split QM/MM at that bond with some probability → multiple candidate QM regions → deduplicate. Avoids exhaustive enumeration.
- **Transfer to our plan:** Use the same idea to define **ORCA fragments**: from each MD frame, extract one or several spherical (or residue-based) fragments, cap with H (or link atoms), run ORCA on each fragment. You then have QM energies/forces **per fragment**, not for the full box. Force-field fitting can target fragment-level residuals or be used to parametrize residue/atom types that appear in many fragments.

### 2.2 Systematic protein partitioning for fragmentation (Wolter et al., *J. Chem. Theory Comput.* 2020)

**Paper:** Wolter, M.; von Looz, M.; Meyerhenke, H.; Jacob, C. R. “Systematic Partitioning of Proteins for Quantum-Chemical Fragmentation Methods Using Graph Algorithms.”  
**arXiv:** 2010.02832

- **Idea:** Protein = **graph** (residues = nodes). Edges are weighted by **estimated fragmentation error** when cutting that bond. Use **graph partitioning** to get a near-optimal partition for a given **max fragment size** and a **local target quantity** (e.g. interaction energy at active site, spectrum of a chromophore).
- **Methods referenced:**  
  - **MFCC** (molecular fractionation with conjugate caps): partition into residue-sized fragments + capping; can add electrostatic embedding.  
  - **FMO** (fragment molecular orbital): overlapping fragments + more sophisticated embedding.
- **Result:** Graph-based partitioning **reduces fragmentation error** compared to naïve fixed-size fragmentation when varying fragment size.
- **Transfer to our plan:** Before running ORCA on many MD frames, define a **single** (or a few) partitioning scheme(s) for your protein: e.g. residue-based graph, edge weights from connectivity or importance to the property you care about. Then for each MD snapshot, **extract the same fragment definitions** (with coordinates from the snapshot), cap, run ORCA. Fitting can use fragment energies/forces consistently across frames; graph-based cuts minimize systematic error for the chosen target (e.g. binding site, specific backbone region).

### 2.3 Fragment-based data for force-field parametrization (ByteFF, Zheng et al. 2024)

**Paper:** Zheng, T. et al. “Data-Driven Parametrization of Molecular Mechanics Force Fields for Expansive Chemical Space Coverage.”  
**arXiv:** 2408.12817 | *Chem. Sci.* 2024

- **Idea:** Don’t run QM on full drug-like molecules everywhere. Use **2.4M optimized molecular fragment geometries** (with Hessians) + **3.2M torsion profiles** at DFT (B3LYP-D3/def2). Train a **GNN** to predict MM parameters (bonded + nonbonded) from structure → **ByteFF** (Amber-compatible).
- **Transfer to our plan:** You don’t have to fit the *full* box at once. You can:
  - Define **chemical fragments** (residue types, side chains, small motifs) that appear in your protein/system.
  - Run **ORCA** on these fragments in many conformations (from MD or from sampling).
  - Feed **fragment** energies/forces (and optionally Hessians/torsion scans) into ForceBalance or a lightweight optimizer to fit **transferable** parameters (e.g. residue parameters, atom types). Full-system accuracy then comes from using these parameters in OpenMM for the full box.

### 2.4 QM force fields and linear scaling (Giese et al., *J. Chem. Theory Comput.*)

**Paper:** Giese, T. J.; et al. “Recent Advances toward a General Purpose Linear-Scaling Quantum Force Field.”  
**Snippet (Semantic Scholar):** QMFFs divide the system into **fragments**, each treated with QM, with **empirically modeled interactions** between fragments → high accuracy + efficiency, good for extensive configurational sampling.

- **Transfer to our plan:** If long term you want “QM-like” accuracy in MD, you could aim for a **fragment-based potential** (each fragment has a small QM or ML model; interactions fitted). For the *immediate* pipeline (ORCA + classical FF fitting), the takeaway is: **fragment-level QM is a valid and established route** to generating reference data that can drive system-level or transferable force fields.

---

## 3. Recommended transfer to your pipeline

### 3.1 Making ORCA feasible

1. **Define fragments once (graph-based or sphere-based)**  
   - Residue graph + max fragment size (e.g. 1–3 residues, or sphere radius 6 Å) so each fragment is &lt;100–150 atoms.  
   - Optionally use Wolter-style edge weights (e.g. importance to binding site or to the property you want to fit).

2. **Per MD frame (or per selected frame in active learning):**  
   - Extract coordinates for each fragment; **cap** cut bonds (H or link atoms).  
   - Run **ORCA** (single-point or gradient) on each fragment only → get E and F per fragment.  
   - Option A: Store **fragment** E/F and fit a force field that has fragment-level (or residue-type) terms.  
   - Option B: Use fragment results to approximate **local** contributions (e.g. active-site energy) and fit parameters that affect that region.

3. **ORCA interface**  
   - Use **OPI** (ORCA Python Interface) to generate inputs and parse outputs.  
   - Write a small adapter: **fragment geometry → ORCA (via OPI) → qdata-like format** (or per-fragment E/F list). Your optimizer (ForceBalance or lightweight) then reads this as “reference data” (possibly with one target per fragment type or per fragment instance).

### 3.2 Making force-field fitting feasible

1. **Fit transferable parameters, not the full box**  
   - Use fragment QM data to optimize **residue/atom-type** parameters (e.g. in OpenMM XML). ForceBalance already supports parameterize-by-element and by type; you can restrict to parameters that appear in your fragments.  
   - Full MD box then uses the **same** XML; no need to run QM on the full box.

2. **Targets**  
   - **Fragment ab initio target:** One (or more) ForceBalance-style targets where “reference” = ORCA fragment E/F and “model” = OpenMM energy/force on the **same** fragment (same topology, same coordinates).  
   - Optionally **torsion profiles** or **optimized geometries** for representative fragments (as in ByteFF) to improve dihedrals and bonded terms.

3. **Active learning loop with fragments**  
   - Run **OpenMM** on full system (with current FF).  
   - **Select** frames (e.g. uncertainty, diversity).  
   - For selected frames: **partition** into fragments → run **ORCA** on each fragment → append to reference set (per fragment type or per fragment).  
   - **Fit** (ForceBalance or lightweight) on accumulated fragment data → update OpenMM XML.  
   - Repeat.  
   - This keeps ORCA cost per iteration bounded by (number of selected frames × fragments per frame × cost per fragment), all with small fragments.

### 3.3 Practical choices

| Decision | Options |
|----------|--------|
| Fragment definition | Residue-based (MFCC-style) vs. sphere-based (Brunken–Reiher) vs. graph-optimized (Wolter). |
| Capping | H atoms vs. link atoms; keep consistent with how you assign atom types in the FF. |
| What to fit | Only parameters that appear in fragments (e.g. residue types in your protein); shared across full system. |
| Reference data format | Per-fragment qdata.txt (or equivalent) so your optimizer can treat each fragment as a small “molecule” target. |

---

## 4. Summary

- **Fragmentation** (sphere-based, residue-based, or graph-optimized) is a standard way to make **QM tractable** for large systems and to generate **reference data for force-field parametrization** without running QM on the full MD box.
- **ORCA** can be run only on **capped fragments** (&lt;150 atoms); use **OPI** to drive ORCA and parse results.
- **Force-field fitting** can use **fragment-level** energies/forces (and optionally Hessians/torsions) to optimize **transferable** parameters (residue/atom types) in OpenMM XML; the full system then uses the same parameters.
- **Active learning** with OpenMM + ORCA + fitting works with this design: select frames → fragment → ORCA on fragments → add to reference set → fit → update FF → repeat. Fragmentation keeps both ORCA and fitting feasible.

---

## 5. References (key papers)

1. Brunken, C.; Reiher, M. *J. Chem. Theory Comput.* **2021**, *17*, 3799–3815. “Automated Construction of Quantum-Classical Hybrid Models.” arXiv:2102.09355.  
2. Wolter, M.; von Looz, M.; Meyerhenke, H.; Jacob, C. R. *J. Chem. Theory Comput.* **2020**, *16*, 6358–6376. “Systematic Partitioning of Proteins for Quantum-Chemical Fragmentation Methods Using Graph Algorithms.” arXiv:2010.02832.  
3. Zheng, T. et al. *Chem. Sci.* **2024**. “Data-Driven Parametrization of Molecular Mechanics Force Fields for Expansive Chemical Space Coverage.” arXiv:2408.12817.  
4. Giese, T. J.; Huang, M.-B.; Chen, H.; York, D. *J. Chem. Theory Comput.* (Special Issue: Beyond QM/MM: Fragment Quantum Mechanical Methods). “Recent Advances toward a General Purpose Linear-Scaling Quantum Force Field.”
