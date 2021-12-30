#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <malloc.h>
#include <stdlib.h>

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

	for(int i = 0; i < b*r*c; i++)
		send_first[i] = drand48();

	for(int i = 0; i < b*r*c; i++)
		send_second[i] = drand48();
    }
    
    MPI_Scatter(send_first,b/size*r*c,MPI_DOUBLE,&localfirst,b/size*r*c,MPI_DOUBLE,0,matrix_comm);
    MPI_Scatter(send_second,b/size*r*c,MPI_DOUBLE,&localsecond,b/size*r*c,MPI_DOUBLE,0,matrix_comm);

    for (int i = 0; i < b/size*r*c; i++)
	sum[i] = localfirst[i] + localsecond[i];

    MPI_Gather(&sum,b/size*r*c,MPI_DOUBLE,recv_buf,b/size*r*c,MPI_DOUBLE,0,matrix_comm);

    MPI_Barrier(matrix_comm);
    end = MPI_Wtime();

    MPI_Finalize( );
    
    if (rank == 0) {
	printf("Runtime, %dD, %d, %d, %d, %f\n", dim, b, r, c, end-start);
    }

    return 0;
}
