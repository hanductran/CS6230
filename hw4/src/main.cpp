#include <iostream>
#include <cstdlib>
#include <random>
#include "quicksort.h"
#include <time.h>
#include <omp.h>

int main() {
    const int NUM = 20; // number of elements of A
    float *A = (float *) malloc(sizeof(float) * NUM);

    for (int i = 0; i < NUM; i++) {
        A[i] = std::rand() % NUM;
    }
    std::cout << "input array initialized" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << A[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;

    double t1= omp_get_wtime();
    quicksort(A, 0, NUM - 1);
    double t2= omp_get_wtime();
    printf("Time taken: %.6fs\n", (t2-t1));

    std::cout << "array sorted" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << A[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;
    delete [] A;
    return 0;
}