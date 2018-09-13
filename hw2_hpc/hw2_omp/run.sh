#!/bin/sh
export OMP_NUM_THREADS=1
./build/hw2_omp
export OMP_NUM_THREADS=2
./build/hw2_omp
export OMP_NUM_THREADS=4
./build/hw2_omp
export OMP_NUM_THREADS=8
./build/hw2_omp
export OMP_NUM_THREADS=16
./build/hw2_omp
export OMP_NUM_THREADS=32
./build/hw2_omp
export OMP_NUM_THREADS=64
./build/hw2_omp
export OMP_NUM_THREADS=128
