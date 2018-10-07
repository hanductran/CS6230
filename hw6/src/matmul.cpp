#include "matmul.h"
#include <omp.h>
#include <iostream>

int matmul(double* A, double* B, double* C, unsigned long int m, unsigned long int k, unsigned long int n) {
    // compute matrix multiplication: C += A*B
    // size of matrices:
    // A=(m,k)
    // B=(k,n)
    // C=(m,n)




    for(unsigned long int i=0;i<m;i++)
        for(unsigned long int j=0;j<n;j++)
        {
            double sum=0;
#pragma omp parallel for default(none) firstprivate(i,j) shared(A,B,k,m,n) reduction(+:sum)
            for(unsigned long int p=0;p<k;p++)
                sum+=A[i*k+p]*B[p*n+j];

            C[i*n+j]+=sum;
        }



    /*#pragma omp parallel
    {
        unsigned int crow, ccol;// arow, acol, brow, bcol, aid, bid;

        unsigned int rank=omp_get_thread_num();;
        const int size = omp_get_num_threads();

        const unsigned long ibegin = ((m*n)*rank)/size;
        const unsigned long iend = ((m*n)*(rank+1))/size;

        //for (unsigned i=0; i<m*n; i++) {
        for (unsigned long i = ibegin; i < iend; i++) {
            crow = (i-(i%n))/n;
            ccol = i%n;
            for (unsigned j=0; j<k; j++){
                C[i] += A[(crow * k + j)]*B[(j * n + ccol)];
            }
        }
    }*/

    return 0;
}
