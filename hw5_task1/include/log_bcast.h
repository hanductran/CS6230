//
// Created by Han Tran on 8/25/18.
//

#ifndef HW5_LOG_BCAST_H
#define HW5_LOG_BCAST_H

#include <mpi.h>
#include <mpi.h>
#include <iostream>
#include <math.h>

int log_bcast(void* data, int count, MPI_Datatype data_type, int root, MPI_Comm comm);

#endif