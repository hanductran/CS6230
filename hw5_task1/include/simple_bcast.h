//
// Created by Han Tran on 8/25/18.
//

#ifndef HW5_SIMPLE_BCAST_H
#define HW5_SIMPLE_BCAST_H

#include <mpi.h>

int simple_bcast(void* data, int count, MPI_Datatype data_type, int root, MPI_Comm comm);

#endif