#include <iostream>
#include <cstdlib>
#include <random>
#include "quicksort.h"

int main() {
    const int NUM = 10; // number of elements of A
    float *A = (float *) malloc(sizeof(float) * NUM);

    for (int i = 0; i < NUM; i++) {
        A[i] = std::rand() % 100;
    }
    std::cout << "input array initialized" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << A[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;
    quicksort(A, 0, NUM - 1);
    std::cout << "array sorted" << std::endl;
    for (unsigned int i = 0; i < NUM; i++) {
        std::cout << A[i] << " ,";
        if (i == 10) std::cout << std::endl;
    }
    std::cout << std::endl;
    return 0;
}