#include <iostream>
#include <math.h>
#include <random>
#include <omp.h>
#include <mpi.h>
//#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]) {

    int rank, size;
    int j, rc;
    int x, xt;
    int partner;
    MPI_Status Stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        double dcheck = log2((double)size);
        int icheck = dcheck;
        if (icheck != dcheck) {
            printf("The number of processes must be a power of 2...\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(0);
        }
    }

    // generate random number for each task
    srand(time(NULL) + rank );
    x = (size-rank)*std::rand() % 100;

    // print initial array
    if (rank == 0)
        printf("Initial array:\n");
    for (int i=0; i<size; i++) {
        if (i == rank) {
            printf("rank= %d, x= %d\n", rank, x);
        }
        // This is just to make the printed data come out in the order of rank for easy check
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // bitonic sort
    int d = log2(size);
    for (int i = 0; i<d; ++i) {
        for (int j = i; j >= 0; --j) {
            partner = rank ^ (1 << j);
            //if (rank.bit(i + 1) == rank.bit(j))
            if ((rank & (1 << (i + 1))) / (1 << (i + 1)) == (rank & (1 << j)) / (1 << j))
            {
                // exchange x with partner, keep the min to replace x
                MPI_Sendrecv(&x, 1, MPI_INT, partner, 1, &xt, 1, MPI_INT, partner, 1, MPI_COMM_WORLD, &Stat);
                if (xt < x) {
                    x = xt;
                }
            } else {
                // exchange x with partner, keep the max to replace x
                MPI_Sendrecv(&x, 1, MPI_INT, partner, 1, &xt, 1, MPI_INT, partner, 1, MPI_COMM_WORLD, &Stat);
                if (xt > x) {
                    x = xt;
                }
            }
        } // i
    } // j

    if (rank == 0)
        printf("Sorted array:\n");
    for (int i=0; i<size; i++) {
        if (i == rank) {
            printf("rank= %d, x= %d\n", rank, x);
        }
        // This is just to make the printed data come out in the order of rank for easy check
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
