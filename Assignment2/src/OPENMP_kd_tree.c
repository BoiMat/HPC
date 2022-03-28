#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#if defined(_OPENMP)
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +     \
                     (double)myts.tv_nsec * 1e-9)

#else

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)
#endif


#if !defined(DOUBLE_PRECISION)
#define float_t float
#else
#define float_t double
#endif

#if !defined(NDIM)
#define NDIM 2
#endif

#define CACHE_LINE 64

typedef unsigned int uint;

typedef struct {
	float_t x[NDIM];
} kpoint;

typedef struct kdnode{
	kpoint split;
	struct kdnode *left, *right;
	uint axis;
} kdnode;

void create_2darray(kpoint *array, uint size);
void swap(kpoint *x, kpoint *y);
kpoint* choose_splitting_point(kpoint *start, kpoint *end, uint axis);	
kdnode* build_kdtree (kpoint *t, uint N, uint axis, uint ndim );
void treeprint(kdnode *root, int level);

int main(int argc, char **argv) {

	struct timespec ts, myts;
	uint threads = 1;
	uint axis = 0;
	uint N = argc >= 2 ? atoi(argv[1]) : 100000000;
	uint c = argc >= 3 ?  atoi(argv[2]) : 0;

	kpoint *data = (kpoint*) malloc (N * sizeof(kpoint));

	create_2darray(data, N);

	#if defined(_OPENMP)
	#pragma omp parallel
	threads = omp_get_num_threads();
	#endif

	double start = CPU_TIME;
	
	struct kdnode *root;

	#pragma omp parallel
	#pragma omp master
	root = build_kdtree(data, N, axis, NDIM);

	double end = CPU_TIME;
	printf("OpenMP: Size = %d\t Threads = %d\t Runtime = %g sec\n", N, threads, end-start);

	if (c != 0)
		treeprint(root, 0);
	
	free(data);

	return 0;
}

void create_2darray(kpoint *array, uint size)
{
	#if defined(_OPENMP)
	#pragma omp parallel
	{
		int me = omp_get_thread_num();
    	short int seed = time(NULL) % ( (1 << sizeof(short int))-1 );
    	short int seeds[3] = {seed-me, seed+me, seed+me*2};

   		#pragma omp for
    	for ( uint i = 0; i < size; i++ ) {
			array[i].x[0] = erand48( seeds );
			array[i].x[1] = erand48( seeds );
		}
	}
	#else 
	{
		long int seed = time(NULL);
		srand48(seed);
		
		printf("seed is % ld\n", seed);
		
		for ( uint i = 0; i < size; i++ ) 
		{
			array[i].x[0] = drand48();
			array[i].x[1] = drand48();
		}
	}    
	#endif
}

void swap(kpoint *x, kpoint *y) 
{
	kpoint t = *x;
	*x = *y;
	*y = t;
}

kpoint* choose_splitting_point(kpoint *start, kpoint *end, uint axis) 
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

kdnode *build_kdtree (kpoint *t, uint N, uint axis, uint ndim ) {
	
	kpoint *n;
	kdnode *node = (kdnode*) malloc(sizeof(kdnode));

	if(!N) return 0;

	if(N==1)
	{
		node->split.x[0] = t->x[0];
		node->split.x[1] = t->x[1];
		node->left = NULL;
		node->right = NULL;
		node->axis = axis;
		return 0;
	}

	n = choose_splitting_point(t, t + N, axis);

	node->split.x[0] = n->x[0];
	node->split.x[1] = n->x[1];
	node->axis = axis;

	axis = (axis + 1) % ndim;
	
	if (N > CACHE_LINE/sizeof(kpoint))
	{
		#pragma omp task
		node->left = build_kdtree(t, n - t, axis, ndim);
		
		#pragma omp task
		node->right = build_kdtree(n + 1, t + N - (n + 1), axis, ndim);
	}
	else
	{
			node->left = build_kdtree(t, n - t, axis, ndim);
			node->right = build_kdtree(n + 1, t + N - (n + 1), axis, ndim);
	}
	return node;
}

void treeprint(kdnode *root, int level)
{
        if (root == NULL)
                return;
        for (int i = 0; i < level; i++)
                printf(i == level - 1 ? "|-" : "  ");
        printf("(%f,%f)\n", root->split.x[0], root->split.x[1]);
        treeprint(root->left, level + 1);
        treeprint(root->right, level + 1);
}