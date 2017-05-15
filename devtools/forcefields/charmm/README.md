This directory contains files and scripts needed to convert the CHARMM forcefield to ffxml files

## Manifest
* ```toppar/``` - All files taken and unzipped from (http://mackerell.umaryland.edu/charmm_ff.html) Jan 2016.
* ```ffxml/``` - charmm36 ffxml
* ```charmm36.yaml``` - yaml file needed for convert_charmm.py script. Specifies files to include and exclude.
* ```convert_charmm.py``` - script to convert charmm top and par files to ffxml.

### Dependencies
* ParmEd

---

Notes on files that were excluded from conversion.

There are two glycolipid stream files with duplicate dihedrals with different values. According to discussion with
Alex MacKerell, the carb glycolipid file should be used so the lipid glycolipid stream file was excluded.

```toppar_all36_prot_aldehydes.str``` and ```toppar_all36_na_modifications.str``` have different values for the angle of
atom types O CD CT2. These files should not be used in the same system so both were excluded. A new atom type is needed
to correct this. If  needed, CGenFF can be used for aldehydes or the user can convert these files at their own risk.


