#rm -rf *.arc *.dyn
#~/bin/tinkerOMM2019/dynamic_omm dimer.xyz -key tinker.key 1 1 0.001 4 298 1 
rm -rf *.arc *.dyn
../dynamic_omm.x dimer.xyz -key tinker.key 1 1 0.001 4 298 1 
