#!/bin/csh -f
# run charmm test cases

#charmm executable
setenv charmm /home/alex/charmm/code/c39b2/exec/gnu/charmm

setenv outdir  .

#mkdir -p $outdir

$charmm mindr:0 mini:0 < test_drude_all_2013c.inp > $outdir/test_drude_all_2013c.out 

#grep for crashed jobs
grep '! I I I I I !' $outdir/*.out
grep '\   XXX   /' $outdir/*.out
grep 'CHECK THOSE INPUTS, ACE' $outdir/*.out

