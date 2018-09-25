#include <iostream>
#include <mpi.h>
#include "simple_bcast.h"
#include "log_bcast.h"
#include <math.h>
#include <time.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    clock_t t1, t2;
    clock_t tmax;


    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int root=0;

    int bcast_data = rank;


    t1 = clock();
    //MPI_Bcast( &bcast_data, 1, MPI_INT, root, MPI_COMM_WORLD);
    //simple_bcast(&bcast_data, 1, MPI_INT, root, MPI_COMM_WORLD);
    log_bcast(&bcast_data, 1, MPI_INT, root, MPI_COMM_WORLD);
    t2 = clock() - t1;

    //printf ("p= %d, MPI_Bcast took %d clicks (%f seconds).\n",rank, t2,((float)t2)/CLOCKS_PER_SEC);

    MPI_Reduce((unsigned long *)&t2, (unsigned long *)&tmax, 1, MPI_UNSIGNED_LONG, MPI_MAX, root, MPI_COMM_WORLD);
    //MPI_Reduce((long double *)&t2, (long double *)&tmax, 1, MPI_LONG_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

    if (rank==root) {
        std::cout<<"Max time = "<< ((float)tmax)/CLOCKS_PER_SEC <<std::endl;
    }




    //std::cout<<"rank: "<<rank<<" bcast_data: "<<bcast_data<<std::endl;

    MPI_Finalize();

    return 0;

}