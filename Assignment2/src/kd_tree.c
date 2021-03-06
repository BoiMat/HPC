#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#if !defined(DOUBLE_PRECISION)
#define float_t float
#else
#define float_t double
#endif

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)

#if !defined(NDIM)
#define NDIM 2
#endif

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

int main(int argc, char **argv) {

	struct timespec ts, myts;
	uint threads = 1;
	uint axis = 0;
	uint N = argc >=2 ? atoi(argv[1]) : 100000000;

	kpoint *data = (kpoint*) malloc (N * sizeof(kpoint));

	create_2darray(data, N);

	double start = CPU_TIME;
	
	struct kdnode *root;

	root = build_kdtree(data, N, axis, NDIM);

	double end = CPU_TIME;
	printf("Serial: Size = %d\t Threads = %d\t Runtime = %f sec\n", N, threads, end-start);
	
	free(data);

	return 0;
}

void create_2darray(kpoint *array, uint size)
{
	long int seed = time(NULL);
	srand48(seed);
		
	//printf("seed is % ld\n", seed);
		
	for ( int i = 0; i < size; i++ ) 
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

	node->left = build_kdtree(t, n - t, axis, ndim);
	node->right = build_kdtree(n + 1, t + N - (n + 1), axis, ndim);

	return node;
}
