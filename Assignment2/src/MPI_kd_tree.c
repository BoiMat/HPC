#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#if !defined(DOUBLE_PRECISION)
#define float_t float
#define MPI_FLOAT_T MPI_FLOAT
#else
#define float_t double
#define MPI_FLOAT_T MPI_DOUBLE
#endif

#if !defined(NDIM)
#define NDIM 2
#endif

typedef struct {
	float_t x[NDIM];
} kpoint;

typedef struct kdnode{
	kpoint split;
	struct kdnode *left, *right;
	uint axis;
} kdnode;

int rank, size;
void create_2darray(kpoint *array, int size);   				// function to generate the array with random numbers in [0,1] using the function drand48()
void swap(kpoint *x, kpoint *y);						// function to swap two kpoints 
kpoint *choose_splitting_point(kpoint *start, kpoint *end, int axis);		// function that performs a quickselect according to the current axis and returns the median
kdnode *build_kdtree(kpoint *t, int N, int axis, int ndim, int depth);		// function that builds the tree recursively

int main(int argc, char **argv) {

	double start, end;
	int last, axis = 0;
	int N = argc >=2 ? atoi(argv[1]) : 100000000;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Datatype p_kpoint;	
	MPI_Type_contiguous(2, MPI_FLOAT_T, &p_kpoint);
	MPI_Type_commit(&p_kpoint);

	int seed = time(NULL);
	srand48(seed);

	if (rank == 0) {							// The master process generates the array and starts building the tree

		kpoint *data = (kpoint *) malloc (N * sizeof(kpoint));

		create_2darray(data, N);
		
		#ifdef DEBUG 
		printf("Here proc %d, the 2D-array is\n", rank);
		for (int k=0; k<N; k++)
			printf("(%f,%f)\n", data[k].x[0], data[k].x[1]); 
		#endif	
		
    		start = MPI_Wtime();

		kdnode *Tree;
		Tree = build_kdtree(data, N, axis, NDIM, -1);
	}
	
	// The other processes wait until they receive a part of the array, and then start constructing the subtrees
	else
	{
		kpoint *array;
		int depth = floor(log2(rank));
		int *info = (int*) malloc(2 * sizeof(int));
		int sendr = rank - pow(2,depth);
		
		#ifdef DEBUG 
		printf("Here proc %d, waiting to receive the subarray from proc %d\n", rank, sendr);
		#endif			

		MPI_Recv(info, 2, MPI_INT, sendr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		array = (kpoint*) malloc (info[0]*sizeof(kpoint));
		MPI_Recv(array, info[0], p_kpoint, sendr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		#ifdef DEBUG 
		printf("Here proc %d, depth is %d. Started building the tree with a subarray of %d points\n", rank, depth, info[0]);
		#endif

		kdnode *subtree;
		subtree = build_kdtree(array, info[0], info[1], NDIM, depth);

		//free(info);
		//free(array);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == 0)
	{
		end = MPI_Wtime();
		printf("MPI: Size = %d\t Procs = %d\t Runtime = %f sec\n", N, size, end-start);
	}
	MPI_Finalize();
	return 0;
}

void create_2darray(kpoint *array, int size)
{
	for(int i = 0; i < size; i++)
	{
		array[i].x[0] = drand48();
		array[i].x[1] = drand48();
	}
}

void swap(kpoint *x, kpoint *y)
{
	kpoint t = *x;
	*x = *y;
	*y = t;
}

kpoint *choose_splitting_point(kpoint *start, kpoint *end, int axis) 
{		
	if (end <= start) 
		return NULL;
	
	if (end == start + 1)
		return start;

	kpoint *p, *store, *median = start + (end - start) / 2;
	float_t pivot;

	while (1) {
		pivot = median->x[axis];

		swap(median, end - 1);
		for (store = p = start; p < end; p++) {
			if (p->x[axis] < pivot) {
				if (p != store) 
					swap(p, store);
				store++;
			}
		}
		swap(store, end - 1);

		if (store->x[axis] == median->x[axis])
			return median;
		
		if (store > median) 
			end = store;
		else 
			start = store;
	}
}

kdnode *build_kdtree(kpoint *t, int N, int axis, int ndim, int depth)
{
	MPI_Datatype p_kpoint;	
	MPI_Type_contiguous(2, MPI_FLOAT_T, &p_kpoint);
	MPI_Type_commit(&p_kpoint);
	
	if(!N) return 0;

	kdnode *node = (kdnode*) malloc(sizeof(kdnode));
	kpoint *n;
	
	// if there is only one point left, the returned node will not have children
	if(N==1) {
		node->split.x[0] = t->x[0];
		node->split.x[1] = t->x[1];
		node->left = NULL;
		node->right = NULL;
		node->axis = axis;
		return 0;
	}

	n = choose_splitting_point(t, t + N, axis);	// find the median

	node->split.x[0] = n->x[0];
	node->split.x[1] = n->x[1];
	node->axis = axis;

	axis = (axis + 1) % ndim;
	depth++;

	int rightdim = t + N - (n + 1);

	if (pow(2,depth) >= size)			// if the number of nodes in the next level (2^depth) is greater than the number of processes,
	{						// recursively build the tree without sending the right subarray anymore. 
		
		#ifdef DEBUG 
		printf("Here proc %d, depth is %d. Building the children nodes with a subarray of %d points\n", rank, depth, recv);
		#endif
		
		node->left = build_kdtree(t, n - t, axis, ndim, depth);

		if(rightdim > 0)
			node->right = build_kdtree(n + 1, rightdim, axis, ndim, depth);
		else
			node->right = NULL;
	}
	else 						// if the number of nodes in the next level (2^depth) is smaller than the number of processes,
	{						// the current process will keep the left subarray and send the right one to a free process.
		int *info = (int*) malloc(2 * sizeof(int));
		int recv = pow(2,depth) + rank;
		info[0] = rightdim; 
		info[1] = axis;
		
		#ifdef DEBUG 
		printf("Here proc %d, depth is %d. Sending to proc %d a subarray with %d points\n", rank, depth, rightdim, recv);
		#endif

		MPI_Send(info, 2, MPI_INT, recv, 0, MPI_COMM_WORLD);
		MPI_Send((n+1), rightdim, p_kpoint, recv, 0, MPI_COMM_WORLD);

		node->left = build_kdtree(t, n - t, axis, ndim, depth);
		
		free(info);
	}
	return node;
}
