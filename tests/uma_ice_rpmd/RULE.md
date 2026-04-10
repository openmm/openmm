# UMA Ice RPMD Output Location

**Never use `/tmp` for test or benchmark outputs.** Outputs in `/tmp` are ephemeral and may be deleted on reboot.

## Required

- Always pass `--output` (or `-o`) with a **persistent path** when running:
  - `test_uma_ice_rpmd.py` → `--output tests/uma_ice_rpmd` or similar
  - `benchmark_ase_reference.py` → `-o tests/uma_ice_rpmd` or similar
- Default `output_dir='.'` writes to the current working directory; ensure CWD is the project root or a known persistent location.

## Example

```bash
cd /media/extradrive/Trajectories/openmm
python tests/uma_ice_rpmd/test_uma_ice_rpmd.py --input tests/uma_ice_rpmd/water1.xyz --output tests/uma_ice_rpmd
python tests/uma_ice_rpmd/benchmark_ase_reference.py -o tests/uma_ice_rpmd
```

Outputs will be in `tests/uma_ice_rpmd/`.
