#!/bin/sh
current_dir=$(pwd)
hdf5_path=${HDF5_HOME%/} #HDF5 library for output

cd ../CODE_DIR/athena
#configure
./configure.py --prob=binary_disk_Zhao --eos=adiabatic --cxx=icc -mpi -hdf5 --hdf5_path=$hdf5_path --mpiccmd=mpiicpc
#./configure.py --prob=binary_disk_ZhaoTEST --eos=adiabatic --cxx=gcc -debug
#compile
make clean
make -j
#copy executable
cp bin/athena $current_dir

