
#include "partition.h"

 unsigned int partition(void* base, size_t num, size_t size, int (*compare)(const void*, const void*)) {

    // allocate memory
    void* pivot=malloc(size);
    void* temp=malloc(size);
    size_t i =0;

    memcpy(pivot,(char*)base + (num-1)*size, size);
    //std::cout<<"pivot: "<<(*(int*)pivot);

    for (unsigned int j = 0; j < (num-1); j++) {
        if(compare((char*)base + j*size, pivot)<0)
        {
            memcpy(temp, (char*)base + i*size, size);
            memcpy((char*)base + i*size, (char*)base + j*size, size);
            memcpy((char*)base + j*size, temp, size);
            i++;
        }
    }
    memcpy(temp, (char*)base + (num-1)*size, size);
    memcpy((char*)base + (num-1)*size, (char*)base + i*size, size);
    memcpy((char*)base + i*size, temp, size);

    free(pivot);
    free(temp);

    return i;
 }