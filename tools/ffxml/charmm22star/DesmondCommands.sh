/opt/schrodinger/utilities/pdbconvert -ipdb ../native.pdb -omae native.mae
mae2dms native.mae native.dms
viparr ./native.dms out.dms -d /opt/schrodinger/desmond-v31023/data/viparr/ff3/charmm22star/ --without-constraints

#Note: you also have to manually fix the box size in out.dms before starting...

desmond --include ../../md.cfg --cfg boot.file=./out.dms
