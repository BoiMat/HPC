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
    int dim = atoi(argv[1]);
    int b = atoi(argv[2]);
    int r = atoi(argv[3]); 
    int c = atoi(argv[4]);

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //if (rank == 0) {
    //	printf("enter the blocks, rows and columns: \n");
    //	scanf("%d %d %d", &b, &r, &c);
    //}

    //MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *dims;
    dims = malloc(dim*sizeof(int));
    for (int i=0; i<dim; i++)
	    dims[i] = 0;
	
    MPI_Dims_create(size, dim, dims);

    int *periods;
    periods = malloc(dim*sizeof(int));
    for (int i=0; i<dim; i++)
	    periods[i] = false;

    int reorder = true;

    MPI_Comm matrix_comm;
    MPI_Cart_create( MPI_COMM_WORLD, dim, dims, periods, reorder, &matrix_comm );

    int my_rank;
    MPI_Comm_rank( matrix_comm, &my_rank );

    MPI_Comm_size( matrix_comm, &size );

    int my_coords[dim];
    MPI_Cart_coords( matrix_comm, my_rank, dim, my_coords );
    
    MPI_Barrier(matrix_comm);
    start = MPI_Wtime();

    //printf("[MPI process, old %d, new %d] I am located at (%d, %d). This is b: %d, this is r: %d, this is c: %d\n", rank, my_rank, my_coords[0],my_coords[1], b,r,c);
    
    srand48(my_rank);
    
    int n;
    double localfirst[b/size*r*c];
    double localsecond[b/size*r*c];
    double sum[b/size*r*c];

    double *send_first, *send_second, *recv_buf;
  
    if (rank==0){
	
	send_first = (double *) malloc(b*r*c*sizeof(double));
	send_second = (double *) malloc(b*r*c*sizeof(double));
	recv_buf = (double *) malloc(b*r*c*sizeof(double));

    	build_matrix(send_first, b*r*c, 0);
    	build_matrix(send_second, b*r*c, 0);
    }
    
    MPI_Scatter(send_first,b/size*r*c,MPI_DOUBLE,&localfirst,b/size*r*c,MPI_DOUBLE,0,matrix_comm);
    MPI_Scatter(send_second,b/size*r*c,MPI_DOUBLE,&localsecond,b/size*r*c,MPI_DOUBLE,0,matrix_comm);
    
    sum_matrix(localfirst, localsecond, sum, b/size*r*c, 0);
    
    MPI_Gather(&sum,b/size*r*c,MPI_DOUBLE,recv_buf,b/size*r*c,MPI_DOUBLE,0,matrix_comm);
    
    MPI_Barrier(matrix_comm);
    end = MPI_Wtime();

    if (rank == 0) {
	    free(send_first);
	    free(send_second);
	    free(recv_buf);
	    printf("Memory deallocated\n");
    }
	
    MPI_Finalize( );
    
    if (rank == 0) {
	printf("Runtime, %dD, %d, %d, %d, %f\n", dim, b, r, c, end-start);
    }

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
