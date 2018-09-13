#include <iostream>
#include <mpi.h>
#include <omp.h>

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Status status;
    MPI_Request request_send;
    MPI_Request request_recv;
    MPI_Request req;

    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &size );

    int val = 2*rank;


    //std::cout << "Process  " << rank << " val =  " << val << std::endl;
    printf("Before sending, process %d val = %d\n", rank, val);

    if (rank == size-1) {
        //MPI_Send(&val, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        //MPI_Isend(&val, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &request_send);
        MPI_Isend(&val, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &req);
    } else {
        //MPI_Send(&val, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD);
        //MPI_Isend(&val, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD,&request_send);
        MPI_Isend(&val, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD,&req);
    }



    if (rank == 0) {
        //MPI_Recv(&val, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, &status);
        //MPI_Irecv(&val, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, &request_recv);
        MPI_Irecv(&val, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, &req);
    } else {
        //MPI_Recv(&val, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);
        //MPI_Irecv(&val, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &request_recv);
        MPI_Irecv(&val, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &req);
    }


    MPI_Wait(&req, &status);

    //std::cout << "Process  " << rank << " val =  " << val << std::endl;
    printf("After receiving, process %d val = %d\n", rank, val);


    MPI_Finalize();

    return 0;
}