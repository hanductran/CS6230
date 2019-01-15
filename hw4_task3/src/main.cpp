#include <iostream>
#include <cstdlib>
#include <random>
#include "Rec_Scan.h"
#include <time.h>
#include <omp.h>

int main() {

    const unsigned int NUM = (1<<22);
    long *A = (long *) malloc(sizeof(long) * NUM);
    long *B = (long *) malloc(sizeof(long) * NUM);

    for (int i = 0; i < NUM; i++) {
        A[i] = std::rand() % 10;
    }

    /*std::cout << "initial array" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << A[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;*/

    double t1= omp_get_wtime();
    B = Rec_Scan(A, NUM);
    double t2= omp_get_wtime();
    printf("Time taken: %.6fs\n", (t2-t1));

    /*std::cout << "array after scan" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << B[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;*/

    delete [] A;
    delete [] B;

    return 0;
}