# CHARMM36 paramter and topology files

These files were taken from (http://mackerell.umaryland.edu/charmm_ff.shtml)

# files that should be excluded from conversion.
# toppar_all36_lipid_list.str just lists the lipid str files (used in charmm to load all stream file)

# There are two glycolipid stream files with duplicate dihedrals with different values. According to discussion with
# Alex MacKerell, the carb glycolipid file should be used.

# These two files have different values for the angle of atom types O CD CT2. These files should not be used in the same
# system. A new atom type is needed to correct this. If needed, CGenFF can be used for aldehydes or the user can convert
# these files at their own risk.
