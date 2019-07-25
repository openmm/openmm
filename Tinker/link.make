#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##             (Intel Fortran for Linux Version)             ##
#  ##                                                           ##
#  ###############################################################
#
#
ifort -O0 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a ; strip analyze.x
ifort -O0 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testgrad.x
ifort -O0 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a ; strip timer.x
