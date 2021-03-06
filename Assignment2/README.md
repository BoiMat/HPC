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

## Run

To run the codes with a 10000 points array, 4 parallel workers:

Serial:

` ./kd_tree.x 10000	` 

MPI:

` mpirun -np 4 mpi_tree.x 10000 /sys/ 2> /dev/null `

OpenMP:

` export OMP_NUM_THREADS=4 `

` ./openmp_tree.x 10000 `

An additional argument can be added to print the final tree, e.g.:

` ./openmp_tree.x 10000 1`

## Output
The output will be in the form:

` Size = 10000     Threads = 12    Runtime = 0.00418282 sec `

## Timings

The timings both in a Strong scaling and Weak scaling scenario are present in the "timings" folder in the form of .csv files.
In the strong scaling scenario they report the number of parallel workers (MPI processes or OpenMP threads), the average timing and related standard deviation for an array of 10^8 points.
In the weak scaling scenario also the size of the array is present.
