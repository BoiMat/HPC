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
