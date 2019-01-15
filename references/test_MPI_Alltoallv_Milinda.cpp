#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#define COMM MPI_COMM_WORLD

#define ARRAY_SIZE 1000

/*
 * Test of MPI_Alltoallv
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
    const int stride = 5;
    const int count = 2;

    int rbuf[count];

    MPI_Init(&argc,&argv);
    MPI_Comm_size(COMM, &numtasks);
    MPI_Comm_rank(COMM, &rank);

    /*sendbuf = (int *)malloc(numtasks*stride*sizeof(int));
    displs = (int *)malloc(numtasks* sizeof(int));
    scounts = (int *)malloc(numtasks*sizeof(int));


    root = 0;
    for (int i=0; i<numtasks; ++i) {
        displs[i] = i*stride;
        scounts[i] = count;
    }

    // initialize data for scattering
    for (int i=0; i<numtasks; ++i) {
        for (int j=0; j<count; j++)
            sendbuf[i*stride + j] = rank;
        for (int j=0; j<(stride-count); j++) {
            sendbuf[i*stride + count + j] = rank+1000;
        }
    }

    MPI_Scatterv(sendbuf, scounts, displs, MPI_INT, rbuf, count, MPI_INT, root, COMM);

    // checking if each process receives correct data
    for (int i=0; i<numtasks; i++) {
        if (rank == i)
            for (int j=0; j<count; j++) {
                printf("Task %d rbuf[%d]= %d\n", rank, j, rbuf[j]);
            }
    }*/


    unsigned int data[ARRAY_SIZE];
    for(unsigned int i=0;i<ARRAY_SIZE;i++)
       data[i]=(unsigned int)rand()%numtasks;

    unsigned int *sendCounts=new unsigned int [numtasks];
    unsigned int *recvCounts=new unsigned int[numtasks];

    unsigned int *sendOffset=new unsigned int[numtasks];
    unsigned int *recvOffset=new unsigned int[numtasks];

    for(unsigned int i=0;i<numtasks;i++)
    {
        sendCounts[i]=0;
        recvCounts[i]=0;

    }

    for(unsigned int i=0;i<ARRAY_SIZE;i++) {
        for (unsigned int p = 0; p < numtasks; p++) {
            if (data[i] <= p)
                sendCounts[p]++;
        }
    }

    sendOffset[0]=0;
    recvOffset[0]=0;
    for(unsigned int i=1;i<numtasks;i++)
        sendOffset[i]=sendOffset[i-1]+sendCounts[i];


    std::vector<unsigned int > sendBuffer;
    sendBuffer.resize(sendOffset[numtasks-1]+sendCounts[numtasks-1]);

    for(unsigned int i=0;i<numtasks;i++)
        sendCounts[i]=0;

    for(unsigned int i=0;i<ARRAY_SIZE;i++) {
        for (unsigned int p = 0; p < numtasks; p++) {
            if (data[i] <= p)
            {
                sendBuffer[sendCounts[p]]=data[i];
                sendCounts[p]++;
            }

        }
    }


    MPI_Alltoall(sendCounts,1,MPI_INT,recvCounts,1,MPI_INT,COMM);

    for(unsigned int i=1;i<numtasks;i++)
        recvOffset[i]=recvOffset[i-1]+recvCounts[i];


    std::vector<unsigned int > recvBuffer;
    recvBuffer.resize(recvCounts[numtasks-1]+recvOffset[numtasks-1]);

    MPI_Alltoallv(&(*(sendBuffer.begin())),(int *)sendCounts,(int *)sendOffset,MPI_INT,&(*(recvBuffer.begin())),(int *)recvCounts,(int *)recvOffset,MPI_INT,COMM);


    for(unsigned int i=0;i<recvBuffer.size();i++)
    {
        std::cout<<"rank: "<<rank<<" value : "<<recvBuffer[i]<<std::endl;
    }

    delete [] sendCounts;
    delete [] recvCounts;
    delete [] sendOffset;
    delete [] recvOffset;



    /*free(sendbuf);
    free(displs);
    free(scounts);*/




    MPI_Finalize();
}
