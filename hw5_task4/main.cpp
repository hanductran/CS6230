
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define COMM MPI_COMM_WORLD
/*
 * Test of MPI_Scatterv
 * root sends out sets of "count" integers from sendbuf[] to each process (each process receives different set)
 * but the sets are "stride" apart in the sendbuf (in this test, strike = 12, count = 10).
 * For testing, the sendbuff is initilized as sendbuf = [[set 1] [set 2] ... [set p]] where each set contains 12 elements
 * but the only first 10 elements are identical (which equals to the rank of process that it is going to receive) while the last 2 elements are equal to 2*rank.
 * The print out to check if each process receives exactly the first 10 elements from each set.
 *
 * This test refers to the example shown at: https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node72.html
 */
int main (int argc, char *argv[])
{
    int numtasks, rank, i;
    int *sendbuf;
    int root, *displs, *scounts;
    const int stride = 12;
    const int count = 10;

    int rbuf[count];

    MPI_Init(&argc,&argv);
    MPI_Comm_size(COMM, &numtasks);
    MPI_Comm_rank(COMM, &rank);

    sendbuf = (int *)malloc(numtasks*stride*sizeof(int));
    displs = (int *)malloc(numtasks* sizeof(int));
    scounts = (int *)malloc(numtasks*sizeof(int));


    root = 0;
    for (int i=0; i<numtasks; ++i) {
        displs[i] = i*stride;
        scounts[i] = count;
    }

    // initialize data in root for scattering
    if (rank==root)
        for (int i=0; i<numtasks; ++i) {
            for (int j=0; j<count; j++)
                sendbuf[i*stride + j] = i;
            for (int j=0; j<(stride-count); j++) {
                sendbuf[i*stride + count + j] = 2*i;
            }
        }

    MPI_Scatterv(sendbuf, scounts, displs, MPI_INT, rbuf, count, MPI_INT, root, COMM);

    // checking if each process receives correct data
    for (int i=0; i<numtasks; i++) {
        if (rank == i)
            for (int j=0; j<count; j++) {
                printf("Task %d rbuf[%d]= %d\n", rank, j, rbuf[j]);
            }
    }
    free(sendbuf);
    free(displs);
    free(scounts);
    MPI_Finalize();
}