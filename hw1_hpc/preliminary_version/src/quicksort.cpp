#include "quicksort.h"
#include "partition.h"


void quicksort(float * A, int lo, int hi) {
    unsigned int p;

    if (lo < hi) {
        p = partition(A, lo, hi);
        quicksort(A, lo, p - 1);
        quicksort(A, p + 1, hi);
    }
}