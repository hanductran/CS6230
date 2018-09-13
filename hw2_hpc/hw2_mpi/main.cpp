#include <iostream>
#include <mpi.h>
#include <omp.h>

int main(int argc, char* argv[]) {
    // argc is always 1, argv[0] is the name of the program, argv[1]...argv[n] are arguments could be passed to the program
    std::cout<<"argc:"<<argc<<" argv: "<<argv[0]<<std::endl;
    MPI_Init(&argc, &argv);

    int rank_mpi, size_mpi;

    MPI_Comm_rank ( MPI_COMM_WORLD, &rank_mpi );
    MPI_Comm_size ( MPI_COMM_WORLD, &size_mpi );

//    std::cout << "Hello from process  " << rank_mpi << " of " << size_mpi << std::endl;

    #pragma omp parallel
    {
        int size_omp, rank_omp;
        rank_omp = omp_get_thread_num();
        size_omp = omp_get_num_threads();
//        std::cout << "Hello from " << rank_omp
//                    << "of" << size_omp << "processes."
//                    << std::endl;
        printf("Hello World from (rank_mpi,rank_omp)= %d %d \
               of (size_mpi,size_omp)= %d %d\n",rank_mpi, rank_omp, size_mpi, size_omp);
    };

    MPI_Finalize();

    return 0;
}