#include "matmul.h"
#include <omp.h>
#include <iostream>

int matmul(double* A, double* B, double* C, int m, int k, int n) {
    // compute matrix multiplication: C += A*B
    // size of matrices:
    // A=(m,k)
    // B=(k,n)
    // C=(m,n)
    unsigned crow, ccol, arow, acol, brow, bcol, aid, bid;

#pragma omp parallel
    {
        int rank;
        const int size = omp_get_num_threads();
        rank = omp_get_thread_num();
        const unsigned long ibegin = ((m*n)*rank)/size;
        const unsigned long iend = ((m*n)*(rank+1))/size;

        //for (unsigned i=0; i<m*n; i++) {
        for (unsigned long i = ibegin; i < iend; i++) {
            crow = (i-(i%n))/n;
            ccol = i%n;
            for (unsigned j=0; j<k; j++){
                arow = crow;
                acol = j;
                aid = arow * k + acol;
                brow = j;
                bcol = ccol;
                bid = brow * n + bcol;
                C[i] += A[aid]*B[bid];
            }
        }
    };

    return 0;
}
