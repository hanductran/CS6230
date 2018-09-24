#include "log_bcast.h"


int log_bcast(void* data, int count, MPI_Datatype data_type, int root, MPI_Comm comm) {

    int rank, size;
    int srank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Status status;

    // ceiling of log(p)
    unsigned n;
    if (size == 1) {
        n = 0;
    } else {
        n = (unsigned)(ceil(log2(size)));
    };
    //std::cout<<"n= "<<n<<std::endl;

    for (unsigned i=0; i<n; i++) {
        // shift rank to the left by root elements so that root starts as "zero" rank
        if (rank < root) {
            srank = rank + (size - root);
        } else {
            srank = rank - root;
        }

        if (srank < (1<<i)) {
            // send
            if ((srank + (1<<i)) < size) {
                // destination exists, then send
                if ((srank + (1<<i)) < (size-root)) {
                    //printf("1.sending: i= %d, rank= %d, srank= %d, send to %d\n", i, rank, srank, srank + (1<<i) + root);
                    MPI_Send(data, count, data_type, srank + (1<<i) + root, 0, comm);
                } else {
                    //printf("2.sending: i= %d, rank= %d, srank= %d, send to %d\n", i, rank, srank, srank + (1<<i) - size + root);
                    MPI_Send(data, count, data_type, srank + (1<<i) - size + root, 0, comm); // 0 is the tag
                }
            }
        } else if (srank < (1<<(i+1))){
            //receive
            if ((srank - (1<<i)) < (size-root)) {
                //printf("3.receiving: i= %d, rank= %d, srank= %d, receive from %d\n", i, rank, srank, srank - (1<<i) + root);
                MPI_Recv(data, count, data_type, srank - (1<<i) + root, 0, comm, &status);
            } else {
                //printf("4.receiving: i= %d, rank= %d, srank= %d, receive from %d\n", i, rank, srank, srank - (1<<i) - size + root);
                MPI_Recv(data, count, data_type, srank - (1<<i) - size + root, 0, comm, &status);
            }

        }
    }

    return 0;
}
