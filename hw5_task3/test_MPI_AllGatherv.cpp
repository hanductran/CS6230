#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define COMM MPI_COMM_WORLD
/*
 * Test of MPI_AllGatherv
 * Every process sends its set of 2 integers (value = rank) to rbuf[] of all processes.
 * The sets are received at locations "stride" apart (in this test, strike = 5).
 *
 * The print out to check if the root receives correct data at correct locations
 */
int main (int argc, char *argv[])
{
    int numtasks, rank, i;

    int root, *rbuf, *displs, *counts;
    const int stride = 5;
    const int count = 2;

    int sendbuf[count];


    MPI_Init(&argc,&argv);
    MPI_Comm_size(COMM, &numtasks);
    MPI_Comm_rank(COMM, &rank);

    rbuf = (int *)malloc(numtasks*stride*sizeof(int));
    displs = (int *)malloc(numtasks* sizeof(int));
    counts = (int *)malloc(numtasks*sizeof(int));

    root = 0;

    for (int i=0; i<numtasks; ++i) {
        displs[i] = i*stride;
        counts[i] = count;
    }

    // Initialize rbuf at each process
    for (int i=0; i<numtasks*stride; i++)
        rbuf[i] = 10000;

    // Initialize data at each process for gathering
    for (int i=0; i<count; ++i) {
        sendbuf[i] = rank;
    }

    MPI_Allgatherv(sendbuf, count, MPI_INT, rbuf, counts, displs, MPI_INT, COMM);

    // checking if the root receive correct data
    for (int i=0; i<numtasks; i++) {
        for (int j=0; j<stride; j++)
            printf("Task %d rbuf[%d]= %d\n", rank, i*stride+j , rbuf[i*stride+j]);
    }

    free(rbuf);
    free(counts);
    free(displs);

    MPI_Finalize();

}
