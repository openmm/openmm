## OpenMM Examples

This directory contains various simple examples demonstrating the use of the
OpenMM Python and C++ APIs.

### Python API examples

The examples in `python-examples` demonstrate how to get a simple simulation up
and running  using the OpenMM Python API application layer.  The way of using
OpenMM illustrated here is primarily geared towards running biomolecular
simulations, but can be used for any kind of simulation that is to be set up by
reading from:
- Amber format files (see `python-examples/simulateAmber.py`)
- CHARMM format files (see `python-examples/simulateCharmm.py`)
- GROMACS format files (see `python-examples/simulateGromacs.py`)
- PDB (Protein Data Bank) files (to be used with OpenMM-compatible force fields; see `python-examples/simulatePdb.py`)

These examples can also be found in the
[Running Simulations chapter](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)
of the user guide, along with explanations of how they work.

### C++ API examples

The examples in `cpp-examples` demonstrate the use of OpenMM's C++ API, and also
show how to use its C and Fortran bindings.  For more information, see
[the README file](cpp-examples/README.md) in this subdirectory.

### Extras

You can also find:
- [Extra utility scripts](extras/README.md) in `extras`
- A suite of benchmarks for OpenMM in `benchmarks`

### Cavity-Coupled MD Examples

The `cavity/` directory contains examples of cavity-coupled molecular dynamics simulations:
- `cavity/dimer_system/` - Simple two-component dimer system with cavity coupling
- `cavity/water_system/` - Flexible TIP4P-Ew water with cavity coupling
- `cavity/protein_system/` - Protein system data and examples

These examples demonstrate:
- Cavity-molecule coupling effects
- IR spectrum analysis with MESA
- Rabi splitting in molecular systems
- Performance benchmarking

See the README files in each subdirectory for detailed usage instructions.

### ML Potential Examples

The `ml/` directory contains examples using machine learning potentials:
- `ml/uma_ice_rpmd/` - UMA model examples for ice RPMD simulations

These examples demonstrate integration of ML potentials with OpenMM for enhanced accuracy in molecular dynamics simulations.
