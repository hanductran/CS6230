#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
//#include <vector>
#include <iostream>

#define COMM MPI_COMM_WORLD

/*
 * Test of MPI_Alltoallv
 * Each process sends to process rank with rank elements starting from index rank of the sbuf.
 * Each process will receive exactly rank elements from all process, and put to the position so that
 * after receiving from all processes, the rbuf will get (rank*nprocs) which are filled from the starting index of rbuf.
 *
 * This test refers to the example shown at: http://mpi.deino.net/mpi_functions/MPI_Alltoallv.html
 */
 
int main (int argc, char *argv[])
{
    int numtasks, rank, *p;

    int *sbuf, *rbuf;
    int *scounts, *rcounts;
    int *sdispls, *rdispls;

    MPI_Init(&argc,&argv);
    
    MPI_Comm_size(COMM, &numtasks);
    MPI_Comm_rank(COMM, &rank);

    sbuf = (int *)malloc(numtasks*numtasks*sizeof(int));
    rbuf = (int *)malloc(numtasks*numtasks*sizeof(int));

    scounts = (int *)malloc(numtasks*sizeof(int));
    rcounts = (int *)malloc(numtasks*sizeof(int));

    sdispls = (int *)malloc(numtasks*sizeof(int));
    rdispls = (int *)malloc(numtasks*sizeof(int));

    // initialize data to send
    for (int i=0; i<numtasks*numtasks; ++i) {
        sbuf[i] = i;
        rbuf[i] = -1000000;
    }

    for (int i=0; i<numtasks; ++i) {
        scounts[i] = i;
    }
    // rcounts which is the transpose of scounts
    for (int i=0; i<numtasks; ++i) {
        rcounts[i] = rank;
    }
    for (int i=0; i<numtasks; ++i) {
        sdispls[i] = i*(i+1)/2;
    }
    for (int i=0; i<numtasks; ++i) {
        rdispls[i] = rank*i;
    }

    MPI_Alltoallv(sbuf, scounts, sdispls, MPI_INT, rbuf, rcounts, rdispls, MPI_INT, COMM);

    // checking if each process receives correct data
    for (int i=0; i<numtasks*numtasks; i++) {
            printf("rank= %d, rbuf[%d]= %d\n",rank, i,rbuf[i]);
    }
    free(sbuf);
    free(scounts);
    free(sdispls);
    free(rbuf);
    free(rcounts);
    free(rdispls);

    MPI_Finalize();
    return 0;
}
