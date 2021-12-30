# Assigment 1 FHPC course 2021/2022


## Section  1: MPI programming

Here all the commands to run the programs are listed. It is assumed that the programmer is inside an ORFEO computational node while running this commands.

### Ring

First thing is to load the essential modules for MPI programming:

` module load openmpi-4.1.1+gnu-9.3.0 `

Then the codes can be compiled.

Odd - even processes approach:

` mpicc ring_Odd-Even.c -o ring_Odd-Even.x `

Non-blocking approaches:

` mpicc ring_Send-Recv.c -o ring_Send-Recv.x `

` mpicc ring_Isend-Irecv.c -o ring_Isend-Irecv.x `

` mpicc ring_Isend-Recv.c -o ring_Isend-Recv.x `

Then they can be run (here with 4 processes):

` mpirun -np 4 ./ring_Odd-Even.x `

` mpirun -np 4 ./ring_Send-Recv.x `

` mpirun -np 4 ./ ring_Isend-Irecv.x `

` mpirun -np 4 ./ring_Isend-Recv.x `

### Matrix

The Matrix code is compiled the same way:

` mpicc sum3Dmatrix.c -o sum3Dmatrix.x`

The code takes as arguments of the main function three integer representing:

- Number of dimension extensions for the virtual topology
 
- The three dimensions for the matrices

So the command would look something like this for a 3D topology and 2400x100x100 matrices:

` mpirun -np 4 ./sum3Dmatrix.x 3 2400 100 100 `

For the other cases present in the GitHub repository (sum3Dmatrix_noVT.c and sum3Dmatrix_serial.c) given the fact that there is not virtual topology the corresponding argument in the command is not needed and the code can be run with:

` mpirun -np 4 ./sum3Dmatrix_noVT.x 2400 100 100 `
