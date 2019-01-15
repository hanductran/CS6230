//
// Created by Han Tran on 8/25/18.
//

#ifndef HW4_HPC_PREFIXSUM_H

#include <iostream>
#include <omp.h>
#include <math.h>

//long* Rec_Scan(const long* A, unsigned int n);
template <class T>
T* Rec_Scan(const T* A, unsigned int n);

#endif

//long* Rec_Scan(const long * A, unsigned int n) {
template <class T>
T* Rec_Scan(const T* A, unsigned int n) {

    //long* S=new long[n];
    T* S=new long[n];

    const long d=log2(n);
    unsigned idx;

    // upsweep
    #pragma omp parallel
    {
        #pragma omp for
        for (unsigned i = 0; i < n; i++)
            S[i] = A[i];
    }
    for (unsigned l=1; l < d+1; l++) {
        #pragma omp parallel
        {
            #pragma omp for
            for (unsigned i=0; i < (1<<(d-l)); i++) {
                idx = (1<<l)*(i+1) - 1;
                S[idx] = S[idx] + S[idx - (1<<(l-1))];
            };
        }

    };

    //downsweep
    for (unsigned l=d-1; l > 0; l--) {
        #pragma omp parallel
        {
            #pragma omp for
            for (unsigned i=0; i < (1<<(d-l)) -1; i++) {
                idx = (1<<l)*(i+1) - 1;
                S[idx + (1<<(l-1))] += S[idx];
            };
        }

    }
    return S;

}