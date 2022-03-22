#!/bin/bash

cd src

gcc kd_tree.c -o kd_tree.x -lm
mpicc MPI_kd_tree.c -o mpi_tree.x -lm
gcc -o openmp_tree.x -fopenmp OPENMP_kd_tree.c

N=1000
P=4

export OMP_NUM_THREADS=$P

./kd_tree.x $N	
mpirun -np $P mpi_tree.x $N /sys/ 2> /dev/null
./openmp_tree.x $N