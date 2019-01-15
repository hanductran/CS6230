#include <iostream>
#include <cstdlib>
#include <random>
#include "quicksort.h"
#include <time.h>
#include <omp.h>

int sort_test(const long* base, int num) {
    unsigned int i;
    for (unsigned int i = 0; i < (num-1); i++) {
        if (base[i] > base[i+1]) {
            return 1;
        }
    }
    return 0;
}

int main() {
    const int NUM = 4000000; // number of elements of A
    long *A = (long *) malloc(sizeof(long) * NUM);

    for (int i = 0; i < NUM; i++) {
        A[i] = std::rand() % 1000;
    }
    /*std::cout << "input array initialized" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << A[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;*/

    double t1= omp_get_wtime();
    quicksort(A, 0, NUM - 1);
    double t2= omp_get_wtime();
    printf("Time taken: %.6fs\n", (t2-t1));

    /*std::cout << "array sorted" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << A[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;*/

    // test if the array is sorted
    int success = sort_test(A, NUM);
    if (success == 0) {
        std::cout<<"Succeed! Array is sorted.";
    } else {
        std::cout<<"Fail! Array is not sorted.";
    }
    std::cout<<std::endl;

    delete [] A;
    return 0;
}