#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    double start, end;
    int rank, value, size, itag_r, itag_l, msgleft, msgright, tmpleft, tmpright, false=0, true=1;
    int right_p, left_p;
    int n_iter = 1;
    MPI_Comm   ring_comm;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    MPI_Cart_create( MPI_COMM_WORLD, 1, &size, &true, 1, &ring_comm );
    MPI_Cart_shift( ring_comm, 0, 1, &left_p, &right_p );
    MPI_Comm_rank( ring_comm, &rank );
    MPI_Comm_size( ring_comm, &size );

    MPI_Barrier(ring_comm);
    start = MPI_Wtime();

    MPI_Status stats[4];
    MPI_Request reqs[4];
    
    for (size_t i=0; i<n_iter; i++) {
	    
	    msgleft = rank;
	    msgright = -rank;
	    
	    itag_r = rank*10;
	    itag_l = rank*10;
	    
	    do {
		    MPI_Irecv( &tmpleft, 1, MPI_INT, right_p, MPI_ANY_TAG, ring_comm, &reqs[0] );
		    MPI_Irecv( &tmpright, 1, MPI_INT, left_p, MPI_ANY_TAG, ring_comm, &reqs[1] );
		    
		    MPI_Isend( &msgleft, 1, MPI_INT, left_p, itag_l, ring_comm, &reqs[2] );
		    MPI_Isend( &msgright, 1, MPI_INT, right_p, itag_r, ring_comm, &reqs[3] );
		    
		    MPI_Waitall(4, reqs, stats);
		    
		    itag_l = stats[0].MPI_TAG;
		    itag_r = stats[1].MPI_TAG;
		    
		    msgleft = tmpleft - rank;
		    msgright = tmpright + rank;
	    
	    } while(rank*10 != itag_r);
    }
	
    MPI_Barrier(ring_comm);
    end = MPI_Wtime();
	
    for(int i=0; i<size; i++){
	    if (rank == i)
		    printf("I am processor %d and i have received %d messages. My final messages have tag %d and value %d, %d\n", rank, size*2, itag_r, tmpleft, tmpright);
    }
	
    MPI_Finalize( );

    if (rank == 0)
	    printf("Runtime, %d, %f\n", size, end-start);

    return 0;
}
