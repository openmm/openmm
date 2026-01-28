# Submodule layout for muhammadhasyim

**fairchem** uses your repo: `https://github.com/muhammadhasyim/fairchem.git` — no change needed.

**cace** uses your fork: `https://github.com/muhammadhasyim/cace.git`. LES-BEC content (ferro_PbTiO3, superionic_water, water) has been merged into this repo at top level, so there is no separate LES-BEC submodule.

## Layout

- **cace**: Your fork at `muhammadhasyim/cace`. Top-level dirs include `ferro_PbTiO3/`, `superionic_water/`, `water/` (CACE-LR water RPMD, fits, docs) plus the main CACE package and examples.
- **fairchem**: Your repo at `muhammadhasyim/fairchem` (branch `les_branch`).

## Keeping cace in sync

From the **openmm** repo root, after pulling or changing cace:

```bash
git submodule update --init cace
```

To push new cace commits to your fork:

```bash
cd cace
git remote add myfork https://github.com/muhammadhasyim/cace.git 2>/dev/null || true
git push myfork main
cd ..
```
