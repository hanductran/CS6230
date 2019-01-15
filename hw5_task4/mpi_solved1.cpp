/******************************************************************************
* FILE: mpi_bug1.c
* DESCRIPTION:
*   This program has a bug that causes it to hang.
* AUTHOR: Blaise Barney
* LAST REVISED: 04/13/05
******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

/*
 * The bug comes from the tag in send and receive: rank 0 send with tag = 0 but rank 1 try to receive from 0 with tag = 1.
 * Thus rank 1 will wait forever at MPI_Recv, and consequently will never send to task 0.
 * To fix: make the tag of send and receive to be the same, then it works.
 */

int main (int argc, char *argv[])
{
    int numtasks, rank, dest, tag, source, rc, count;
    char inmsg, outmsg='x';
    MPI_Status Stat;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Task %d starting...\n",rank);

    if (rank == 0) {
        if (numtasks > 2)
            printf("Numtasks=%d. Only 2 needed. Ignoring extra...\n",numtasks);
        dest = rank + 1;
        source = dest;
        tag = rank;

        rc = MPI_Send(&outmsg, 1, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
        printf("Sent to task %d...\n",dest);

        rc = MPI_Recv(&inmsg, 1, MPI_CHAR, source, 1, MPI_COMM_WORLD, &Stat);
        printf("Received from task %d...\n",source);
    }

    else if (rank == 1) {
        dest = rank - 1;
        source = dest;
        tag = rank;

        rc = MPI_Recv(&inmsg, 1, MPI_CHAR, source, 0, MPI_COMM_WORLD, &Stat);
        printf("Received from task %d...\n",source);

        rc = MPI_Send(&outmsg, 1, MPI_CHAR, dest, 1, MPI_COMM_WORLD);
        printf("Sent to task %d...\n",dest);

    }


    if (rank < 2) {
        rc = MPI_Get_count(&Stat, MPI_CHAR, &count);
        printf("Task %d: Received %d char(s) from task %d with tag %d \n",
               rank, count, Stat.MPI_SOURCE, Stat.MPI_TAG);
    }
    MPI_Finalize();
    return 0;
}