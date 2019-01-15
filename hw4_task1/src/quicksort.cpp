#include "quicksort.h"
#include "partition.h"
#include <omp.h>

void quicksort(long * A, int lo, int hi) {
    unsigned int p;

    if (lo < hi) {
        p = partition(A, lo, hi);

        if(p>0)
        {
            #pragma omp task
            quicksort(A, lo, p-1);
        }

        #pragma omp task
        quicksort(A, p+1, hi);

    };
}