#include <iostream>
#include <mpi.h>
#include "simple_bcast.h"
#include "log_bcast.h"
#include <math.h>

int main(int argc, char* argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int root=3;

    int bcast_data=rank;

    MPI_Bcast( &bcast_data, 1, MPI_INT, root, MPI_COMM_WORLD);
    //simple_bcast(&bcast_data, 1, MPI_INT, root, MPI_COMM_WORLD);
    //log_bcast(&bcast_data, 1, MPI_INT, root, MPI_COMM_WORLD);

    std::cout<<"rank: "<<rank<<" bcast_data: "<<bcast_data<<std::endl;


    MPI_Finalize();

    return 0;

}