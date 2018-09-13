#include "quicksort.h"
#include "partition.h"

void quicksort(void* base, size_t num, size_t size, int (*compare)(const void*, const void*)) {
    size_t  p;

    if (num>1){
        p = partition(base, num, size, compare);
        quicksort(base, p, size, compare);
        quicksort((char*)base + (p+1)*size, num-(p+1), size, compare);
    }
}