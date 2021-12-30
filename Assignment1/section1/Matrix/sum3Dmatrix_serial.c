#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <malloc.h>
#include <stdlib.h>

void build_matrix (double *matrix, int size, int print);
void sum_matrix(double *m1, double *m2, double *sum, int size, int print);

int main(int argc, char *argv[])
{
    double start, end;
    int rank, size;

    int b = atoi(argv[1]);
    int r = atoi(argv[2]); 
    int c = atoi(argv[3]);

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    //printf("[MPI process, old %d, new %d] I am located at (%d, %d). This is b: %d, this is r: %d, this is c: %d\n", rank, my_rank, my_coords[0],my_coords[1], b,r,c);
    
    srand48(rank);
    
    double *sum; 
    sum = (double *) malloc(b*r*c*sizeof(double));

    double *send_first, *send_second; //*recv_buf;
  
    if (rank==0){
	
	send_first = (double *) malloc(b*r*c*sizeof(double));
	send_second = (double *) malloc(b*r*c*sizeof(double));

    	build_matrix(send_first, b*r*c, 0);
    	build_matrix(send_second, b*r*c, 0);

        sum_matrix(send_first, send_second, sum, b*r*c, 0);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    if (rank == 0) {
	    free(send_first);
	    free(send_second);
    }
	
    MPI_Finalize( );
    
    if (rank == 0)
	printf("Runtime, %d, %d, %d, %f\n", b, r, c, end-start);
    

    return 0;
}


void build_matrix (double *matrix, int size, int print) {
	int i;
	for(i = 0; i < size; i++) {
		matrix[i] = drand48();
		if (print == 1)
			printf("%f ", matrix[i]);
	}
}


void sum_matrix(double *m1, double *m2, double *sum, int size, int print) {
	int i;

	for (i = 0; i < size; i++) {
		sum[i] = m1[i] + m2[i];
		if (print == 1)
			printf("%f ", sum[i]);
	}
}
