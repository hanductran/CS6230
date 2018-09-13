#!/bin/sh
export OMP_NUM_THREADS=8
mpirun -np 4 ./build/hw2_mpi
