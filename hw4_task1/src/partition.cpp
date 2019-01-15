#include <omp.h>
#include "partition.h"

unsigned int partition(long* A, unsigned int lo, unsigned int hi) {

    long pivot=A[hi];
    unsigned int i=lo;
    float tempo;

    //#pragma omp parallel for default(none) firstprivate(pivot) shared(lo,hi,A,i) private(tempo)
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
