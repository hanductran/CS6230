#include <iostream>
#include <vector>
#include <functional>
#include <omp.h>
#include <cmath>
#include <time.h>
#include "dense_matvec.h"

int main() {

    const unsigned int n = 300000, m=25000;

    double *A = (double *)malloc(sizeof(double)*n*m);
    for (unsigned int i=0; i<n*m; i++) {
        A[i] = (double)std::rand()/(double)(RAND_MAX/5.0);
    }

    double *X = (double *)malloc(sizeof(double)*m);
    for (unsigned int i=0; i<m; i++) {
        X[i] = (double)std::rand()/(double)(RAND_MAX/5.0);
    }

    double *Y = (double *)malloc(sizeof(double)*n);

    double t1= omp_get_wtime();

    dense_matvec(A, X, Y, n, m);

    double t2=omp_get_wtime();
    printf("Time taken: %.6fs\n", (t2-t1));

    delete [] A;
    delete [] X;
    delete [] Y;
    return 0;
}