//
// Created by Han Tran on 9/10/18.
//
#include <omp.h>
#include <iostream>
#include "SpMv.h"
#include <time.h>

int SpMv(double* val, unsigned int* col, unsigned int* rpt, double* x, double* y, unsigned int n) {
    // This functions implements multiplication of a sparse matrix with a vector x, output is vector y
    // Matrix is square of size n by n, which is stored in compressed row storage format using 3 vectors:
    // val[nnz] contains the value of non-zero elements
    // col[nnz] contains the column-index of the non-zero elements
    // rpt[n+1] contains the row-index range of the non-zero elements
    int p = omp_get_max_threads();
    int chunksize = (n/p)?16:1;

    // initialize results
#pragma omp parallel for
    for (unsigned k = 0; k < n; k++) {
        y[k] = 0;
    }

    // matrix-vector multiplication
#pragma omp parallel for schedule(dynamic, chunksize)
    for (unsigned i = 0; i < n; i++) {
        for (unsigned k = rpt[i]; k < rpt[i+1]; k++) {
            y[i] = y[i] + val[k]*x[col[k]];
    }
    }
    return 0;
}
