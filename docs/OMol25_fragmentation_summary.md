# What OMol25 Says About Fragmentation (and How It Compares)

**Short answers:**  
- OMol25 uses **pocket/shell extraction + capping** to get DFT-sized systems from large biomolecules and electrolytes. It does **not** present a general QM fragmentation *method* (like FMO or graph-optimized partitioning).  
- It is **not** claiming a fragmentation technique that is superior to methods in the fragmentation literature (e.g. Wolter et al., Brunken–Reiher). It describes a **practical, scalable protocol** for building a large ML training dataset.

---

## 1. What OMol25 Actually Does (“Fragmentation” in the Paper)

### 1.1 Main idea

OMol25 needs **DFT-feasible system sizes** (up to 350 atoms). Large PDB/MD systems are made tractable by **extracting local regions** and **capping** cut boundaries. The word “fragment” in the paper refers to these **extracted, capped subsystems**, not to a QM fragmentation scheme (e.g. FMO, MFCC).

### 1.2 Biomolecules (sections/omol_dataset.tex, appendix/bio_methods.tex)

- **Protein–ligand (pocket extraction)**  
  - Take the **immediate protein pocket** around the ligand: BioLiP2 receptor residues + any residue with an atom within **2.5 Å** of the ligand.  
  - **Cap** protein: ACE/NMA at N/C termini (or glycine where a single residue was skipped).  
  - Add H, fix Lewis structures, assign charge/spin; keep systems **≤ 350 atoms**.  
  - Run **MD with Cα restraints** so the extracted pocket doesn’t fall apart.  
  - Subsample frames (e.g. first and last of 10) for DFT.

- **“Fragmented protein pockets + ligand” (bio_methods.tex §2.2)**  
  - After docking drug-like molecules into pockets, **further reduce size**: keep **only the two receptor fragments (chains) closest to the ligand**, with at most **5 residues total**.  
  - Use only the **last MD frame**.  
  - Goal: **smaller** systems so more DFT calculations fit the same compute budget.  
  - This is the only place the paper uses “fragmented” in a distinctive way: **distance-based retention** of the 2 closest pocket fragments.

- **Protein–protein core**  
  - Pick buried residues (low SASA); extract all residues with side-chain heavy atoms within **4–4.5 Å**. Cap, then same MD + subsampling.

- **Protein–protein interface**  
  - DIPS-Plus interface residues; extract within **4 Å** (same chain) or **6 Å** (other chain). Cap, MD, subsample.

- **Nucleic acids**  
  - **Motif-based** extraction around a “core” base (several neighborhood types by connectivity/distance). Cap (e.g. O3′ OH); restrain backbone in MD; subsample frames.

### 1.3 Electrolytes (omol_dataset.tex, electro_methods.tex)

- **Shell extraction:** For each ion (or solvent molecule), take all molecules that have **any atom within 3 Å or 5 Å** of the central ion.  
- Discard clusters &gt; 350 atoms.  
- Same idea: **local extraction by radius** → DFT-sized clusters.

### 1.4 Summary of OMol25 “fragmentation”

| Context | What OMol25 does |
|--------|-------------------|
| Protein–ligand | Extract pocket (2.5 Å around ligand) + cap (ACE/NMA or Gly). |
| “Fragmented” pockets | After docking: keep only 2 closest receptor fragments to ligand (≤5 residues). |
| Protein core/interface | Extract by 4–6 Å from selected residues; cap. |
| Nucleic acids | Motif-based neighborhoods; cap. |
| Electrolytes | 3 Å / 5 Å shells around ions/solvent. |

So across the board: **spatial extraction by distance/radius + capping + (for biomolecules) restrained MD** to get stable, DFT-sized inputs. No optimization of *where* to cut to minimize QM error; no FMO/MFCC-style fragmentation algebra.

---

## 2. Does OMol25 Outline a *Superior* Fragmentation Technique?

**No.** The paper does not claim, and does not provide, a fragmentation technique that is superior to established methods.

- **OMol25’s goal:** Build a **large, diverse DFT dataset** for training MLIPs. Extraction + capping is chosen to be **simple, automatable, and scalable** (fixed cutoffs, standard caps, MD with restraints).  
- **Fragmentation literature (e.g. Wolter et al., Brunken–Reiher):** Aims at **accuracy of a QM quantity** (e.g. interaction energy, spectrum) by choosing *how* to partition the system (graph-based minimization of fragmentation error, sphere-based QM regions with rule-based cuts, etc.).  

OMol25 does **not**:

- Optimize cut positions to minimize fragmentation error.  
- Use graph partitioning or FMO/MFCC-style fragmentation.  
- Compare to other fragmentation schemes or claim superiority.  
- Derive or justify cutoffs (2.5 Å, 4 Å, 3/5 Å) from fragmentation theory; they are practical choices for dataset construction.

So: **OMol25’s approach is a good, practical protocol for dataset generation**, not a new or “superior” fragmentation *method* in the sense of the QM fragmentation literature.

---

## 3. How OMol25’s Approach Relates to Your Plan (ORCA + Force-Field Fitting)

- **Reusable ideas:**  
  - **Extract local regions** (pocket, shell, or residue neighborhood) so ORCA only sees &lt;100–350 atoms.  
  - **Cap** (ACE/NMA, H, or link atoms) so fragments are chemically valid.  
  - Optionally **restrain** in MD (e.g. Cα) so fragments don’t fall apart before QM.  
  - You can use **distance/radius rules** (like 2.5 Å, 4 Å, 3/5 Å) as a first pass; for proteins, the “two closest fragments” idea (retain only a few residues near the ligand) is a simple way to get smaller ORCA inputs.

- **What to add from the fragmentation literature for *accuracy*:**  
  - If you care about **minimizing QM error** for a specific quantity (e.g. binding site energy), consider **graph-based partitioning** (Wolter et al.) or **sphere-based QM regions** (Brunken–Reiher) instead of fixed cutoffs only.  
  - OMol25 does not replace that; it gives you a **scalable extraction recipe** you can combine with more refined fragmentation if needed.

---

## 4. References (in the doc you already have)

- **OMol25:** arXiv:2505.08762 (Open Molecules 2025).  
- **Fragmentation methods:** See `docs/Fragmentation_for_ORCA_and_FF_fitting.md` (Brunken–Reiher, Wolter et al., ByteFF, etc.).
