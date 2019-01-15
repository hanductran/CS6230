#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define COMM MPI_COMM_WORLD
/*
 * Test of MPI_Gatherv
 * Every process sends its set of "count" integers to rbuf[] of root.
 * The sets are received at root at locations "stride" apart (in this test, strike = 12).
 * For this test: each process sends 10 of integers with value of its rank to the root. These sets of 10 integers are positioned
 * at [[0...9], [stride ...(stride + 9)], [2*stride ...(2*stride + 9)]...]
 *
 * The print out to check if the root receives correct data at correct locations. It prints out the first element of set
 * and the element right after the end element of the set
 *
 * This test is based on the example shown at: https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node70.html#Node70
 */
int main (int argc, char *argv[])
{
    int numtasks, rank, i;

    int root, *rbuf, *displs, *counts;
    const int stride = 15;
    const int count = 10;

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

    // Initialize rbuf at root
    if (rank==root)
        for (int i=0; i<numtasks*stride; i++)
            rbuf[i] = 10000;

    // Initialize data at each process for gathering
    for (int i=0; i<count; ++i) {
        sendbuf[i] = rank;
    }

    MPI_Gatherv (sendbuf, 10, MPI_INT, rbuf, counts, displs, MPI_INT, root, COMM);

    // checking if the root receive correct data
    if (rank==root)
        for (int i=0; i<numtasks; i++) {
            printf("Task %d rbuf[%d]= %d; rbuf[%d]= %d\n", rank, i*stride , rbuf[i*stride], i*stride + count , rbuf[i*stride + count]);
        }
    free(rbuf);
    free(counts);
    free(displs);

    MPI_Finalize();

}
