#include "simple_bcast.h"
#include <mpi.h>
#include <iostream>

int simple_bcast(void* data, int count, MPI_Datatype data_type, int root, MPI_Comm comm) {

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Status status;
    if (rank == root) {
        // root --> send data to others
        for (int r = 0; r < size; r++)
            if (r != rank)
                MPI_Send(data, count, data_type, r, 0, comm); // 0 is the tag
    } else {
        MPI_Recv(data, count, data_type, root, 0, comm, &status);
    }

    //printf("Rank %d, data = %d\n", rank, (int *)data);
    
    return 0;
}
