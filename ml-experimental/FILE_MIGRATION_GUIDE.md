# File Migration Summary - ml-experimental Reorganization

This document tracks all files moved during the ml-experimental reorganization. Scripts that reference these files may need path updates.

## Moved Files Reference

### Cavity Particle Files

#### Examples (from tests/ to ml-experimental/cavity/examples/)
- `tests/dimer_system/` → `ml-experimental/cavity/examples/dimer_system/`
- `tests/water_system/` → `ml-experimental/cavity/examples/water_system/`
- `tests/protein_system/` → `ml-experimental/cavity/examples/protein_system/`

**Scripts that may reference these:**
- Any scripts importing from `tests.dimer_system` → update to `ml_experimental.cavity.examples.dimer_system`
- Analysis scripts that reference `tests/dimer_system/*.npz` files

#### Utilities (NEW - created)
- **NEW:** `ml-experimental/cavity/utils/cavity_helpers.py` - Common cavity functions extracted
- **NEW:** `ml-experimental/cavity/utils/__init__.py`

**Scripts should now import:**
```python
from ml_experimental.cavity.utils import (
    add_cavity_particle,
    setup_cavity_coupling,
    wavenumber_to_hartree,
)
```

### RPMD Files

#### Documentation (from root to ml-experimental/rpmd/docs/)
- `HYBRID_RPMD_IMPLEMENTATION_COMPLETE.md` → `ml-experimental/rpmd/docs/`
- `HYBRID_RPMD_QUICKSTART.md` → `ml-experimental/rpmd/docs/`
- `HYBRID_RPMD_STATUS.md` → `ml-experimental/rpmd/docs/`
- `HYBRID_RPMD_V2.md` → `ml-experimental/rpmd/docs/`
- `GPU_OPTIMIZED_BUSSI.md` → `ml-experimental/rpmd/docs/`
- `plugins/rpmd/HYBRID_RPMD.md` → `ml-experimental/rpmd/docs/`

#### Examples (from root to ml-experimental/rpmd/examples/)
- `test_hybrid_rpmd_compile.py` → `ml-experimental/rpmd/examples/`
- `test_hybrid_rpmd_gpu.py` → `ml-experimental/rpmd/examples/`
- `test_hybrid_rpmd_simple.py` → `ml-experimental/rpmd/examples/`
- `test_hybrid_rpmd_v2.py` → `ml-experimental/rpmd/examples/`
- `test_hybrid_rpmd_water.py` → `ml-experimental/rpmd/examples/`
- `test_minimal_rpmd.py` → `ml-experimental/rpmd/examples/`
- `test_minimal_rpmd_reference.py` → `ml-experimental/rpmd/examples/`
- `test_rpmd_pile_g.py` → `ml-experimental/rpmd/examples/`

**To run these scripts now:**
```bash
cd ml-experimental/rpmd/examples
python test_hybrid_rpmd_simple.py
```

### ML Potential Files

#### Documentation (from root to ml-experimental/ml/docs/)
- `AVAILABLE_UMA_MODELS.md` → `ml-experimental/ml/docs/`
- `BATCHED_INFERENCE_EXPLAINED.md` → `ml-experimental/ml/docs/`
- `ESEN_BATCHED_RPMD_FIX.md` → `ml-experimental/ml/docs/`
- `IMPLEMENTATION_COMPLETE.md` → `ml-experimental/ml/docs/`
- `INSTALLATION_STATUS.md` → `ml-experimental/ml/docs/`
- `INSTALLATION_TESTING.md` → `ml-experimental/ml/docs/`
- `OMOL25_LES_INTEGRATION_GUIDE.md` → `ml-experimental/ml/docs/`
- `UMA_IMPLEMENTATION_SUMMARY.md` → `ml-experimental/ml/docs/`
- `UMA_INTEGRATION_FINAL_REPORT.md` → `ml-experimental/ml/docs/`
- `UMA_INTEGRATION_LIMITATION.md` → `ml-experimental/ml/docs/`
- `UMA_PYTHONFORCE_SUCCESS.md` → `ml-experimental/ml/docs/`
- `UMA_QUICKSTART.md` → `ml-experimental/ml/docs/`
- `UMA_RPMD_PERFORMANCE_ANALYSIS.md` → `ml-experimental/ml/docs/`

#### Examples (from root/tests to ml-experimental/ml/examples/)
- `test_rpmd_uma.py` → `ml-experimental/ml/examples/`
- `test_uma_pythonforce.py` → `ml-experimental/ml/examples/`
- `test_uma_pythonforce_simple.py` → `ml-experimental/ml/examples/`
- `test_uma_rpmd_pythonforce.py` → `ml-experimental/ml/examples/`
- `test_uma_torch_export.py` → `ml-experimental/ml/examples/`
- `tests/uma_ice_rpmd/` → `ml-experimental/ml/examples/uma_ice_rpmd/`
- `tests/cace-lr_water/` → `ml-experimental/ml/examples/cace-lr_water/`

**To run these scripts now:**
```bash
cd ml-experimental/ml/examples
python test_uma_pythonforce_simple.py

cd ml-experimental/ml/examples/uma_ice_rpmd
python test_uma_ice_rpmd.py
```

### Other Files Moved
- `cursor_specific_particle_rpmdintegrator.md` → `ml-experimental/`
- `cursor_terminal_log_event.md` → `ml-experimental/`

## Files That Remain in Place

### Core Implementation (unchanged)
- `openmmapi/src/CavityForce.cpp`
- `openmmapi/src/CavityParticleDisplacer.cpp`
- `plugins/rpmd/openmmapi/src/RPMDIntegrator.cpp`
- `plugins/rpmd/platforms/common/src/CommonRpmdKernels.cpp`
- `plugins/uma/` (if exists)
- All platform-specific implementations

### Build System (unchanged)
- `CMakeLists.txt` (root)
- All plugin `CMakeLists.txt` files
- Build configuration files

### Tests (C++ unit tests remain)
- `tests/Test*.h` (C++ test headers)
- `tests/Test*.cpp` (C++ test implementations)
- These are part of the build system and must stay

### Submodules (unchanged)
- `fairchem/` submodule
- `cace/` submodule
- `LES-BEC/` submodule

## Import Path Changes Needed

### Python Scripts

**OLD imports (no longer work):**
```python
# These paths no longer exist
from tests.dimer_system import run_simulation
from tests.water_system import analyze_spectrum
```

**NEW imports (use these):**
```python
# Import from ml-experimental structure
from ml_experimental.cavity.examples.dimer_system import run_simulation
from ml_experimental.cavity.examples.water_system import analyze_spectrum

# Or use the utilities
from ml_experimental.cavity.utils import add_cavity_particle
```

### Running Scripts

**OLD way (no longer works):**
```bash
cd /path/to/openmm
python test_hybrid_rpmd_simple.py  # File not found
```

**NEW way (use these):**
```bash
cd /path/to/openmm/ml-experimental/rpmd/examples
python test_hybrid_rpmd_simple.py

# Or from root with module syntax
cd /path/to/openmm
python -m ml_experimental.rpmd.examples.test_hybrid_rpmd_simple
```

### Documentation Links

Documentation files that referenced moved files should update links:

**OLD links:**
```markdown
See [tests/dimer_system/README.md](tests/dimer_system/README.md)
```

**NEW links:**
```markdown
See [dimer_system README](ml-experimental/cavity/examples/dimer_system/README.md)
```

## Scripts Likely Needing Updates

### Analysis Scripts
Any analysis scripts that reference:
- `tests/dimer_system/*.npz`
- `tests/water_system/*.npz`
- `tests/protein_system/*.pdb`

Should update paths to:
- `ml-experimental/cavity/examples/dimer_system/*.npz`
- `ml-experimental/cavity/examples/water_system/*.npz`
- `ml-experimental/cavity/examples/protein_system/*.pdb`

### Build Scripts
Scripts that run tests from root directory should update:
```bash
# OLD
python test_hybrid_rpmd_simple.py

# NEW
python ml-experimental/rpmd/examples/test_hybrid_rpmd_simple.py
```

### CI/CD Pipelines
If automated tests reference these files, update:
- Test file paths
- Import statements in test runners
- Documentation build paths

## Backwards Compatibility

### Symlinks (Optional)
If needed for temporary backwards compatibility, create symlinks:

```bash
cd /path/to/openmm
ln -s ml-experimental/cavity/examples/dimer_system tests/dimer_system
ln -s ml-experimental/rpmd/examples/test_hybrid_rpmd_simple.py test_hybrid_rpmd_simple.py
```

**Note:** Symlinks are a temporary solution. Update actual paths for long-term maintenance.

### Python Package
To enable imports from `ml_experimental`, ensure the root is in PYTHONPATH:

```bash
export PYTHONPATH="/path/to/openmm:$PYTHONPATH"
```

Or add `__init__.py` files:
```bash
touch ml-experimental/__init__.py
touch ml-experimental/cavity/__init__.py
# etc.
```

## Verification Checklist

After reorganization, verify:

- [ ] All example scripts run from their new locations
- [ ] Documentation links resolve correctly
- [ ] Import statements work (test with Python)
- [ ] Build system still compiles (C++ tests)
- [ ] Git history preserved (files copied, not moved with git mv)
- [ ] No broken symlinks
- [ ] CI/CD pipelines updated

## Git Commands for Future Updates

If you need to move additional files with history:

```bash
# Copy with history (requires git >= 2.32)
git mv old/path new/path
git commit -m "Reorganize: move file to new location"

# Or for bulk moves
git mv tests/dimer_system ml-experimental/cavity/examples/
git commit -m "Reorganize: move dimer_system to cavity examples"
```

## Questions?

See:
- [ml-experimental/README.md](../ml-experimental/README.md) - Overview
- [docs/ml-experimental/ARCHITECTURE.md](docs/ml-experimental/ARCHITECTURE.md) - Architecture
- Individual feature documentation in `ml-experimental/{feature}/docs/`
