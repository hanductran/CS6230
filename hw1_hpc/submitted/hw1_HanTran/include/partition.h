//
// Created by Han Tran on 8/25/18.
//

#ifndef HW1_HPC_PARTITION_H
#define HW1_HPC_PARTITION_H

#include <iostream>
#include <cstring>

unsigned int partition(void* base, size_t num, size_t size, int (*compare)(const void*, const void*));

#endif //HW1_HPC_PARTITION_H
