#include <iostream>
#include <mpi.h>
#include "simple_bcast.h"
#include "log_bcast.h"
#include <math.h>
#include <time.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    double t1, t2;
    //clock_t tmax;
    double tmax;


    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int root=0;

    double bcast_data = (double)std::rand()/(double)(RAND_MAX/5.0);


    t1 = MPI_Wtime();//clock();
    log_bcast(&bcast_data, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    //MPI_Bcast( &bcast_data, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    //simple_bcast( &bcast_data, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    t2 = MPI_Wtime() - t1;

    MPI_Reduce((double *)&t2, (double *)&tmax, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

    if (rank==root) {
        printf ("log_bcast: Max time taken: %.8f\n",tmax);
    }

    MPI_Finalize();

    return 0;

}
