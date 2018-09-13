
#include "partition.h"

unsigned int partition(float* A, unsigned int lo, unsigned int hi) {
    float pivot=A[hi], tempo;
    unsigned int i=lo;
    for (unsigned int j = lo; j < hi; j++) {
        if(A[j] < pivot)
        {
            tempo = A[i];
            A[i] = A[j];
            A[j] = tempo;
            i++;
        }
    }
    tempo = A[hi];
    A[hi] = A[i];
    A[i] = tempo;
    return i;
}
