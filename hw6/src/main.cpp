#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include "matmul.h"

int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Comm comm;
    int dim[2], period[2], reorder;
    int coord[2];
    int source, dest;
    int rc;
    double t1, t2, tmax;


    const unsigned int m = atoi(argv[1]); //number of rows of C and A
    const unsigned int n = atoi(argv[2]); //number of columns of C and B
    const unsigned int k = atoi(argv[3]); //number of columns of A and rows of B

    double* A = (double *) malloc(sizeof(double*) * m * k);
    double* B = (double *) malloc(sizeof(double*) * k * n);
    double* C = (double *) malloc(sizeof(double*) * m * n);

    MPI_Status Stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        double dcheck = sqrt((double)size);
        int icheck = dcheck;
        if (icheck != dcheck) {
            printf("The number of processes must be a perfect square...\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(0);
        }
    }

    int q = (int)sqrt((double)size);

    dim[0] = q;
    dim[1] = q;
    period[0] = 1;
    period[1] = 1;
    reorder = 0;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

    // generate matrices A and B for each rank
    for (unsigned int i = 0; i < m*k; i++) {
        A[i] = (rank+1)*(double)std::rand()/(double)(RAND_MAX/5.0);
        //A[i] = (rank+1)*rand() % 10;
        //A[i] = (rank+1)*(i+1);
    }
    for (unsigned int i = 0; i < k*n; i++) {
        B[i] = (rank+1)*(double)std::rand()/(double)(RAND_MAX/5.0);
        //B[i] = (rank+1)*rand() % 10;
        //B[i] = (rank+2)*(2*i+1);
    }

    // initialize C
    for (unsigned int i = 0; i < m*n; i++) {
        C[i] = 0;
    }

    /*printf("rank= %d\n", rank);
    for (unsigned int i=0; i<m; i++) {
        for (unsigned j=0; j<k; j++) {
            printf("A[%d,%d]=%f\n",i,j,A[i*k+j]);
        }
    }
    MPI_Barrier(comm);
    printf("rank= %d\n", rank);
    for (unsigned int i=0; i<k; i++) {
        for (unsigned j=0; j<n; j++) {
            printf("B[%d,%d]=%f\n",i,j,B[i*n+j]);
        }
    }
    MPI_Barrier(comm);*/

    // TIMING:
    t1 = MPI_Wtime();


    //SKEWING:
    // get Cartesian coordinates of rank
    MPI_Cart_coords(comm, rank, 2, coord);

    // shift rows of A by coord[1] elements (dest will be left of rank, source will be right of rank)
    MPI_Cart_shift(comm, 1, coord[0], &dest, &source);
    // send A to dest and receive A from source
    if (coord[0] > 0) {
        double * recvBufferA=new double[m*k];
        MPI_Sendrecv(A, m*k, MPI_DOUBLE, dest, 1, recvBufferA, m*k, MPI_DOUBLE, source, 1, comm, &Stat);
        std::swap(A,recvBufferA);
        delete [] recvBufferA;
    }

    // shift columns of B by coord[1] (dest will be above of rank, source will be below of rank)
    MPI_Cart_shift(comm, 0, coord[1], &dest, &source);
    // send B to dest and receive B from source
    if (coord[1] > 0) {
        double * recvBufferB=new double[k*n];
        MPI_Sendrecv(B, k*n, MPI_DOUBLE, dest, 1, recvBufferB, k*n, MPI_DOUBLE, source, 1, comm, &Stat);
        std::swap(B,recvBufferB);
        delete [] recvBufferB;
    }

    matmul(A, B, C, m, k, n);

    //SHIFTING: q-1 times for matrix multiplication
    for (unsigned int l=1; l<q; l++) {
        // shift rows of A by 1 element (dest will be left of rank, source will be right of rank)
        MPI_Cart_shift(comm, 1, 1, &dest, &source);
        //send A to dest and receive A from source
        double * recvBufferA=new double[m*k];
        MPI_Sendrecv(A, m*k, MPI_DOUBLE, dest, 1, recvBufferA, m*k, MPI_DOUBLE, source, 1, comm, &Stat);
        std::swap(A,recvBufferA);

        // shift columns of B by 1 element (dest will be above of rank, source will be below of rank)
        MPI_Cart_shift(comm, 0, 1, &dest, &source);
        // send B to dest and receive B from source
        double * recvBufferB=new double[k*n];
        MPI_Sendrecv(B, k*n, MPI_DOUBLE, dest, 1, recvBufferB, k*n, MPI_DOUBLE, source, 1, comm, &Stat);
        std::swap(B,recvBufferB);

        delete [] recvBufferA;
        delete [] recvBufferB;

        // multiply and accumulate C += A * B
        matmul(A, B, C, m, k, n);
    }

    // RETURN A to original position
    MPI_Cart_shift(comm, 1, -(coord[0] - 1), &dest, &source);
    double * recvBufferA=new double[m*k];
    MPI_Sendrecv(A, m*k, MPI_DOUBLE, dest, 1, recvBufferA, m*k, MPI_DOUBLE, source, 1, comm, &Stat);
    std::swap(A,recvBufferA);
    delete [] recvBufferA;

    // RETURN B to original position
    MPI_Cart_shift(comm, 0, -(coord[1] - 1), &dest, &source);
    double * recvBufferB=new double[k*n];
    MPI_Sendrecv(B, k*n, MPI_DOUBLE, dest, 1, recvBufferB, k*n, MPI_DOUBLE, source, 1, comm, &Stat);
    std::swap(B,recvBufferB);
    delete [] recvBufferB;

    /*printf("Final result, rank= %d, A=\n", rank);
    for (unsigned int i=0; i<m; i++) {
        for (unsigned j=0; j<k; j++) {
            printf("A[%d,%d]=%f\n",i,j,A[i*k + j]);
        }
    }
    MPI_Barrier(comm);
    printf("Final result, rank= %d, B=\n", rank);
    for (unsigned int i=0; i<k; i++) {
        for (unsigned j=0; j<n; j++) {
            printf("B[%d,%d]=%f\n",i,j,B[i*n + j]);
        }
    }
    MPI_Barrier(comm);
    printf("Final result, rank= %d, C=\n", rank);
    for (unsigned int i=0; i<m; i++) {
        for (unsigned j = 0; j < n; j++) {
            printf("C[%d,%d]=%f\n", i, j, C[i*n + j]);
        }
    }
    MPI_Barrier(comm);*/

    t2 = MPI_Wtime() - t1;
    MPI_Reduce((double *)&t2, (double *)&tmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    if (rank==0) {
        printf("m= %d, k= %d, n= %d, p= %d; max time taken: %.8f\n", m, k, n, size, tmax);
    }

    free(A);
    free(B);
    free(C);
    MPI_Finalize();

    return 0;
}