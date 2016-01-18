# simple script to run water and the different water-ion pairs individually due to a shake error that is being corrected for c35a2.

alias charmm='/raid/lopes/charmm-raid3/c34a2j-pre-release-of-c34b1/c34a2j-fix/exec/gnu/charmm'

charmm count:0 id1:swm4 < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:1 id1:lit < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:2 id1:sod < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:3 id1:pot < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:4 id1:rub < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:5 id1:ces < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:6 id1:flu < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:7 id1:cla < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:8 id1:bro < positive_drude_test_temp.inp > positive_drude_test_temp.out
charmm count:9 id1:iod < positive_drude_test_temp.inp > positive_drude_test_temp.out
