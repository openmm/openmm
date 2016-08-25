**Summary**: call `python amber2omm.py -v`. Requires the `files/` dir. Outputs `ffxml/` containing all converted XMLs and `log.csv` - log of the validation energies.

**General**

Script supports command line arguments, and so `python amber2omm.py -h`:
```
usage: amber2omm.py [-h] [--input INPUT] [--input-format INPUT_FORMAT]
                    [--verbose] [--no-log] [--protein-test] [--nucleic-test]
                    [--protein-ua-test] [--phospho-protein-test] [--gaff-test]

AMBER --> OpenMM forcefield conversion script

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        path of the input file. Default: "files/master.yaml"
  --input-format INPUT_FORMAT, -if INPUT_FORMAT
                        format of the input file: "yaml" or "leaprc". Default:
                        "yaml"
  --verbose, -v         turns verbosity on
  --no-log              turns logging of energies to log.csv off
  --protein-test        validate resulting XML through protein tests
  --nucleic-test        validate resulting XML through nucleic acid tests
  --protein-ua-test     validate resulting XML through united-atom protein
                        tests
  --phospho-protein-test
                        validate resulting XML through phosphorylated protein
                        tests
  --gaff-test           validate resulting XML through small-molecule (GAFF)
                        test
```

With the defaults as set, all you need to do is have the script and the `files/` directory and call `python amber2omm.py`. `-v` or `--no-log` as wanted.

Output:
* `ffxml/` with all output XMLs
* `log.csv` - csv, machine readable log, containing AMBER, OpenMM, and relative (abs(AMBER-OpenMM) / AMBER) energies  - you can see what this looks like here: https://github.com/choderalab/openmm/blob/forcefield_conversion_bonanza/devtools/forcefield-scripts/amber/log.csv . This is overwritten if already present on a new run, use `--no-log` to disable.
* (LeAP is called extensively by the script, and outputs its own `leap.log` too)

`-v` will turn on printing of the progress-tracking comments from the script. Warnings issued by ParmEd (warnings are by default filtered to error, but in some cases this had to be relaxed to 'always' because of the nature of the files converted) are printed always. The output of LeAP is always redirected to `dev/null`.

**YAML input**

By default the script takes a YAML file input, `files/master.yaml` has everything that is going on here.
There's only a few rules to the required structure of the YAML and it will be very easily extendable to future forcefields.

First entry in the YAML must be:
``` yaml
- sourcePackage: AmberTools
  sourcePackageVersion: 15
```

a MODE declaration follows:
``` yaml
- MODE: LEAPRC
```

There are two modes: LEAPRC and RECIPE. LEAPRC - self-explanatory, RECIPE - this is for water-ion conversions, RECIPE because your input files are a mix of `.dat`, `frcmod` and `.lib`, rather than leaprc.

We're in LEAPRC mode - this is used for all protein/nucleic acids forcefield. So all the FFs follow, example:
``` yaml
- Source: leaprc.ff14SB
  Reference:
  - >-
    Maier, J.A., Martinez, C., Kasavajhala, K., Wickstrom, L., Hauser, K.E., and Simmerling, C. (2015).
    ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB. J. Chem. Theory Comput. 11, 3696-3713.
  Test:
  - protein
  - nucleic
```

I think that's self-explanatory, all fields required.

There's an optional `Options` field, which allows changes to parameters if the default is no good for this case:
``` yaml
- Source: leaprc.phosaa10
  Reference:
  - >-
    Steinbrecher, T., Latzer, J., and Case, D.A. (2012). Revised AMBER parameters for bioorganic phosphates.
    J. Chem. Theory Comput. 8, 4405-4412.
  Options:
    filter_warnings: always
    write_unused: True
  Test:
  - protein_phospho
```

And so we do all protein / nucleic acid FFs.

Then we reach water & ions. MODE is changed to RECIPE:
``` yaml
- MODE: RECIPE
```

and extra source package is declared (water models are converted manually - we just supply them as an ffxml in `files/` - I have `tip3p.xml`, `tip4pew.xml`, `spce.xml` from the current distribution of OpenMM with some changes to make them the 'newest' format (e.g. no classes used, only types). These ffxmls are integrated together with the converted ion parameters to make all the output. So we have OpenMM 7.0 listed as the source of these files - hence the extra input).

``` yaml
- sourcePackage2: OpenMM
  sourcePackageVersion2: 7.0
```

Then again, we're going ffXML by ffXML - but we need a few extra fields. Two examples:
A 'standard' file - water model + JC monovalent ions + compromise set +2 ions
``` yaml
- Source:
  - parm/frcmod.ionsjc_tip3p
  - parm/frcmod.ionslrcm_cm_tip3p
  - lib/atomic_ions.lib
  Solvent_source: tip3p.xml
  Solvent: tip3p
  Name: tip3p_standard
  Reference:
  - >-
    Joung, I.S., and Cheatham, Thomas E. (2008).
    Determination of Alkali and Halide Monovalent Ion Parameters for Use in Explicitly Solvated Biomolecular Simulations. J. Phys. Chem. B 112, 9020-9041.
  - >-
    Joung, I.S., and Cheatham, T.E. (2009).
    Molecular Dynamics Simulations of the Dynamic and Energetic Properties of Alkali and Halide Ions Using Water-Model-Specific Ion Parameters. J. Phys. Chem. B 113, 13279-13290.
  - >-
    Li, P., Roberts, B.P., Chakravorty, D.K., and Merz, K.M. (2013).
    Rational Design of Particle Mesh Ewald Compatible Lennard-Jones Parameters for +2 Metal Cations in Explicit Solvent. J. Chem. Theory Comput. 9, 2733-2748.
  - >-
    Jorgensen, W.L., Chandrasekhar, J., Madura, J.D., Impey, R.W., and Klein, M.L. (1983).
    Comparison of simple potential functions for simulating liquid water. The Journal of Chemical Physics 79, 926-935.
  Test:
  - water_ion
```
An 'overloading' set: HFE +2, +3 and +4 ions for tip3p water.

``` yaml
- Source:
  - parm/frcmod.ionslrcm_hfe_tip3p
  - parm/frcmod.ions34lsm_hfe_tip3p
  - lib/atomic_ions.lib
  Standard: tip3p_standard
  Solvent: tip3p
  Name: tip3p_HFE_multivalent
  Reference:
  - >-
    Li, P., Roberts, B.P., Chakravorty, D.K., and Merz, K.M. (2013).
    Rational Design of Particle Mesh Ewald Compatible Lennard-Jones Parameters for +2 Metal Cations in Explicit Solvent. J. Chem. Theory Comput. 9, 2733-2748.
  - >-
    Li, P., Song, L.F., and Merz, K.M. (2015).
    Parameterization of Highly Charged Metal Ions Using the 12-6-4 LJ-Type Nonbonded Model in Explicit Water. J. Phys. Chem. B 119, 883-895.
  Test:
  - water_ion
```

* Source - AMBER input files

* Solvent_source - the water file in `files/` for the standard (i.e. water model containing) XMLs **or** Standard - this is the same as the `Name` field for the appropriate standard (water model containing) XML - we need to know that, because for the 'overloading' sets both that XML and the standard XML need to be loaded for energy testing

* Solvent - this is the name of the solvent, this is necessary to avoid hardcoding of recognition of what solvent you're using from the names of the files etc. - and knowing which solvent you're using is necessary for energy validations.

* Name - the desired name of the ffxml. (For proteins and nucleic this is done by the script, which a product of `leaprc.ff14SB` will call `ff14SB.xml` etc.)

Should you want to provide a different YAML (shorter version, different completely, whatever you want), script will take it with:
```
python amber2omm.py --input name_of_your_yaml.yaml
```

The outputs of any YAML will be written to `ffxml/`.

You can also provide a leaprc of your choosing via:
```
python amber2omm.py --input name_of_your_leaprc --input-format leaprc
```

The output of any leaprc will be written to the `./`.

**A few remarks on the water models and ions**

As I've said before already, the system is:
* water models converted manually, i.e. taken existing files, put them in `files/`, script will merge them with appropriate ions

* we have *standard* sets: `tip3p_standard.xml`, `tip4pew_standard.xml`, `spce_standard.xml`: water model + JC monovalent ions + compromise set +2 ions

* for each water model we then have an HFE and IOD sets - +2/+3/+4 ions, all have templates set to `overload = "1"`. (`tip3p_HFE_multivalent.xml`, `tip3p_IOD_multivalent.xml` etc.)

* usage is to always load in a standard, and then you can overload +2's and add +3 and +4 with the HFE or IOD files

* naming of the water atom types has stayed as was (`tip3p-O`)

* naming of the ion atom types is `name_of_set (dash) amber_atom_type_name`, e.g. `tip3p_standard-Na+`, `tip3p_HFE_multivalent-Zn2+`.
