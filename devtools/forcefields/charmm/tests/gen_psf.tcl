package require psfgen
topology ../toppar/top_all36_cgenff.rtf
topology ../toppar/toppar_water_ions.str

mol load pdb N_N_dimethylaniline.pdb
set DMAN [atomselect top all]

segment DMAN {
pdb N_N_dimethylaniline.pdb
first NONE
last NONE
}
coordpdb N_N_dimethylaniline.pdb DMAN

regenerate dihedrals angles
guesscoord

writepsf N_N_dimethylaniline.psf
writepdb N_N_dimethylaniline_new.pdb
