//
// Created by Han Tran on 9/10/18.
//
#include <omp.h>
#include <iostream>
#include "dense_matvec.h"


int dense_matvec(double* M, double* x, double* y, unsigned int n, unsigned int m) {
    int p = omp_get_max_threads();
    int chunksize = (n/p)?16:1;

#pragma omp parallel for schedule(dynamic, chunksize)
        for (int i=0; i<n; ++i)
//#pragma omp parallel for
            for (int j=0; j<m; ++j)
                y[i] += M[i*m + j] * x[j];
// return value for error codes.
    return 0;
}

