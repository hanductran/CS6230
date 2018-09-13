//
// Created by Han Tran on 8/25/18.
//

#ifndef HW1_HPC_QUICKSORT_H
#define HW1_HPC_QUICKSORT_H

#include <iostream>

void quicksort(void* base, size_t num, size_t size, int (*compare)(const void*, const void*));

#endif //HW1_HPC_QUICKSORT_H
