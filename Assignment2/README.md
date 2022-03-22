# Assigment 2 FHPC course 2021/2022

## Compile

To run the script the GCC compiler is used and both OpenMP and OpenMPI must be installed.
Inside the scripts both the size of the array and the number of parallel workers can be modified.
The single implementations can be compiled with the following commands:

Serial:
` gcc kd_tree.c -o kd_tree.x -lm `

MPI:
` mpicc MPI_kd_tree.c -o mpi_tree.x -lm `

OpenMP:
` gcc -o openmp_tree.x -fopenmp OPENMP_kd_tree.c `

To run the codes:

Serial:
` ./kd_tree.x 10000	`

MPI:
` mpirun -np $P mpi_tree.x $N /sys/ 2> /dev/null `

OpenMP:
` ./openmp_tree.x $N `
