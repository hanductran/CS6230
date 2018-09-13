#include <iostream>
#include <cstdlib>
#include <random>
#include "quicksort.h"

template <typename T>
int compare (const void * a, const void * b)
{
    if (fabs((*(T*)a) - (*(T*)b)) < 1e-20) {
        return 0;
    } else if ((*(T*)a) < (*(T*)b)) {
        return -1;
    } else {
        return 1;
    }
}

int sort_test(const void* base, size_t num, size_t size, int (*compare)(const void*, const void*)) {
    unsigned int i;
    for (unsigned int i = 0; i < (num-1); i++) {
        if (compare((char*)base + (i+1)*size, (char*)base + i*size) < 0) {
            return 1;
        }
    }
    return 0;
}

int main () {
    int i, success;
    float old_sum=0, new_sum=0;
    const int NUM=100;

    // generate array A of random numbers
    float * A = (float * )malloc(sizeof(float)*NUM);
    for (int i = 0; i < NUM; i++) {
        A[i] = (float)std::rand()/(float)(RAND_MAX/5.0);
    }

    //compute the sum of elements of array A (this is for testing later)
    for (int i = 0; i < NUM; i++) {
        old_sum += A[i];
    }

    // sort the input
    quicksort (A, NUM, sizeof(float), compare<float>);

    // print out sorted input
    std::cout<<"Array after sorting:"<<std::endl;
    for (unsigned i=0; i<NUM; i++) {
        std::cout<<A[i]<<" ,";
        if (((i+1) % 10) == 0) std::cout<<std::endl;
    }
    std::cout<<std::endl;

    // test if the array is sorted
    success = sort_test(A, NUM, sizeof(float), compare<float>);
    if (success == 0) {
        std::cout<<"Succeed! Array is sorted.";
    } else {
        std::cout<<"Fail! Array is not sorted.";
    }
    std::cout<<std::endl;

    // test if the sum of sorted array is equal the sum of unsorted array (i.e. the input)
    for (int i = 0; i < NUM; i++) {
        new_sum += A[i];
    }
    if (old_sum == new_sum){
        std::cout<<"Succeed! The 'sum test' is passed."<<std::endl;
    } else {
        std::cout<<"Fail! The 'sum test' is not passed."<<std::endl;
    }

    return 0;
}