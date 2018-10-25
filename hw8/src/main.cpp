#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#define ID(i,j,k)  ((nnode_z+2) * (nnode_y+2) * i + (nnode_z+2) * j + k)


int main(int argc, char *argv[]) {

    int rank, size, p_size;
    MPI_Comm comm;
    MPI_Status Stat;
    int coord[3], period[3];
    int x_left = -1, x_right = -1;     // proc on the left and right in x direction
    int y_left = -1, y_right = -1;     // proc on the left and right in y direction
    int z_left = -1, z_right = -1;     // proc on the left and right in z direction
    int rc;
    double t1, t2, tmax;
    double xmin, xmax;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    int count;
    int ierr;
    int i,j,k;
    double b, x, y, z;

    //const int N = 100;                     // interior nodes in each direction (i.e. total nodes = N + 2)
    const int N = atoi(argv[1]);
    const int L = 1000 ;                      // side length of the (cube) domain
    const double h = double(L)/double(N+1);// element length
    const unsigned int MAX_ITER=1000;
    double* swap;

    double change, my_change;
    double epsilon = 1.0E-03;
    int step = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        double dcheck = cbrt(size);
        p_size = dcheck;
        if (abs(p_size - dcheck) > 0.01) {
            printf("The number of processes must be cube of an integer, program stops.\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(0);
        }
    }

    p_size = cbrt(size);      // number of processes per direction (x, y and z)

    int periods[3] = {0,0,0};
    int dims[3] = {p_size, p_size, p_size};
    //MPI_Dims_create(size, 3, dims);

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm);
    MPI_Cart_coords(comm, rank, 3, coord);

    //printf("rank= %d, coord= %d %d %d\n", rank, coord[0], coord[1], coord[2]);
    //MPI_Barrier(comm);

    double eps = 0.0001;
    double d = ((N +1 + eps) + eps)/(double)(p_size);// average size of each process

    // min & max nodes for x direction =================================
    xmin = (-eps) + (double)(coord[0]*d);
    xmax = xmin + d;
    for (int i = 0; i <= N+1; i++) {
        if (i >= xmin) {
            ixmin = i;
            break;
        }
    }
    for (int i = N+1; i >= 0; i--) {
        if (i < xmax) {
            ixmax = i;
            break;
        }
    }

    // min & max nodes for y direction =================================
    xmin = (-eps) + (double)(coord[1]*d);
    xmax = xmin + d;
    for (int i = 0; i <= N+1; i++) {
        if (i >= xmin) {
            iymin = i;
            break;
        }
    }
    for (int i = N+1; i >= 0; i--) {
        if (i < xmax) {
            iymax = i;
            break;
        }
    }

    // min & max nodes for z direction =================================
    xmin = (-eps) + (double)(coord[2]*d);
    xmax = xmin + d;
    for (int i = 0; i <= N+1; i++) {
        if (i >= xmin) {
            izmin = i;
            break;
        }
    }
    for (int i = N+1; i >= 0; i--) {
        if (i < xmax) {
            izmax = i;
            break;
        }
    }
    //printf("rank= %d, ixmin= %d, ixmax= %d, iymin= %d, iymax= %d, izmin= %d, izmax= %d\n", rank, ixmin, ixmax,
    //        iymin, iymax, izmin, izmax);


    // rank # of left and right in x, y and z directions =================================
    MPI_Cart_shift(comm, 0, 1, &x_left, &x_right);
    MPI_Cart_shift(comm, 1, 1, &y_left, &y_right);
    MPI_Cart_shift(comm, 2, 1, &z_left, &z_right);
    //printf("rank= %d, x_left= %d, x_right= %d\n", rank, x_left, x_right);


    // allocate =================================
    unsigned int nnode_x = ixmax - ixmin +1; // number of nodes in x direction
    unsigned int nnode_y = iymax - iymin +1; // number of nodes in y direction
    unsigned int nnode_z = izmax - izmin +1; // number of nodes in z direction

    //printf("rank= %d, nnode_x= %d, nnode_y= %d, nnode_z= %d\n", rank, nnode_x, nnode_y, nnode_z);
    //printf("rank= %d, xleft= %d, xright= %d, yleft= %d, yright= %d, zleft= %d, zright= %d\n", rank, x_left, x_right,
    //        y_left, y_right, z_left, z_right);
    //printf("rank= %d, ixmin= %d, ixmax= %d, iymin= %d, iymax= %d, izmin= %d, izmax= %d\n", rank, ixmin, ixmax, iymin,
    //        iymax, izmin, izmax);


    // allocate displacement (size = number of nodes + 2 ghost nodes)
    double* u = (double *)malloc(sizeof(double*) * (nnode_x + 2) * (nnode_y + 2) * (nnode_z + 2));
    double* u_new = (double *)malloc(sizeof(double*) * (nnode_x + 2) * (nnode_y + 2) * (nnode_z + 2));
    for (unsigned i = 0; i < (nnode_x + 2) * (nnode_y+2) * (nnode_z+2); i++) {
        u[i] = 0;
        u_new[i] = 0;
    }

    double* sbuff_x_left = (double *)malloc(sizeof(double*) * nnode_y * nnode_z);
    double* sbuff_x_right = (double *)malloc(sizeof(double*) * nnode_y * nnode_z);
    double* rbuff_x_left = (double *)malloc(sizeof(double*) * nnode_y * nnode_z);
    double* rbuff_x_right = (double *)malloc(sizeof(double*) * nnode_y * nnode_z);

    double* sbuff_y_left = (double *)malloc(sizeof(double*) * nnode_x * nnode_z);
    double* sbuff_y_right = (double *)malloc(sizeof(double*) * nnode_x * nnode_z);
    double* rbuff_y_left = (double *)malloc(sizeof(double*) * nnode_x * nnode_z);
    double* rbuff_y_right = (double *)malloc(sizeof(double*) * nnode_x * nnode_z);

    double* sbuff_z_left = (double *)malloc(sizeof(double*) * nnode_x * nnode_y);
    double* sbuff_z_right = (double *)malloc(sizeof(double*) * nnode_x * nnode_y);
    double* rbuff_z_left = (double *)malloc(sizeof(double*) * nnode_x * nnode_y);
    double* rbuff_z_right = (double *)malloc(sizeof(double*) * nnode_x * nnode_y);

    
    //[0,1] for x-left and x-right; [2,3] for y-left and y-right; [4,5] for z-left and z-right
    MPI_Request	send_request[6], recv_request[6];
    MPI_Status status;

    // Timing
    t1 = MPI_Wtime();

    // WHILE (not convergence) LOOP
    do {
    ++step;

    // apply Dirichlet boundary condition =================================
    if (x_left < 0) {
        i = 1;
        for (int j = 1; j <= nnode_y; j++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = 0.0;
            }
        }
    }
    if (x_right < 0) {
        i = nnode_x;
        for (int j = 1; j <= nnode_y; j++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = 0.0;
            }
        }
    }
    if (y_left < 0) {
        j = 1;
        for (int i = 1; i <= nnode_x; i++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = 0.0;
            }
        }
    }
    if (y_right < 0) {
        j = nnode_y;
        for (int i = 1; i <= nnode_x; i++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = 0.0;
            }
        }
    }
    if (z_left < 0) {
        k = 1;
        for (int i = 1; i <= nnode_x; i++) {
            for (int j = 1; j <= nnode_y; j++) {
                u[ID(i,j,k)] = 0.0;
            }
        }
    }
    if (z_right < 0) {
        k = nnode_z;
        for (int i = 1; i <= nnode_x; i++) {
            for (int j = 1; j <= nnode_y; j++) {
                u[ID(i,j,k)] = 0.0;
            }
        }
    }


    /*if (rank == 0) {
        for (int i=1; i <= 1; i++) {
            for (int j = 1; j <= nnode_y; j++) {
                for (int k = 1; k <= nnode_z; k++) {
                    printf("rank= %d, (i,j,k)= %d, %d, %d, u= %f\n", rank, i, j, k, u[ID(i,j,k)]);
                }
            }
        }
    }*/

    // exchange data with the left in x direction=================================
    if (x_left >= 0) {
        // collect data for sending
        count = 0;
        i = 1;
        for (int j = 1; j <= nnode_y; j++) {
            for (int k = 1; k <= nnode_z; k++) {
                sbuff_x_left[count] = u[ID(i,j,k)];
                count = count + 1;
            }
        }
        // send data to x-left proc
        MPI_Isend(sbuff_x_left, nnode_y * nnode_z, MPI_DOUBLE, x_left, 0, comm, &send_request[0]);
        // receive data from x-left proc
        MPI_Irecv(rbuff_x_left, nnode_y * nnode_z, MPI_DOUBLE, x_left, 1, comm, &recv_request[0]);
        //printf("ghost ex\n");
    }

    // exchange data with the right in x direction=================================
    if (x_right >= 0) {
        // collect data for sending
        count = 0;
        i = nnode_x;
        for (int j = 1; j <= nnode_y; j++) {
            for (int k = 1; k <= nnode_z; k++) {
                sbuff_x_right[count] = u[ID(i,j,k)];
                count = count + 1;
            }
        }
        // send data to x-right proc
        MPI_Isend(sbuff_x_right, nnode_y * nnode_z, MPI_DOUBLE, x_right, 1, comm, &send_request[1]);
        // receive data from x-right proc
        MPI_Irecv(rbuff_x_right, nnode_y * nnode_z, MPI_DOUBLE, x_right, 0, comm, &recv_request[1]);
    }


    // exchange data with the left in y direction =================================
    //CHECKING
    /*if (rank == 22) {
        for (unsigned i = ixmin; i <= ixmax; i++) {
            for (unsigned j = iymin; j <= iymax; j++) {
                for (unsigned k = izmin; k <= izmax; k++) {
                    u[index(i,j,k)] = index(i,j,k);
                }
            }
        }
    }*/
    if (y_left >= 0) {
        // collect data for sending
        count = 0;
        j = 1;
        for (int i = 1; i <= nnode_x; i++) {
            for (int k = 1; k <= nnode_z; k++) {
                sbuff_y_left[count] = u[ID(i,j,k)];
                count = count + 1;
            }
        }
        // send to y-left proc
        MPI_Isend(sbuff_y_left, nnode_x * nnode_z, MPI_DOUBLE, y_left, 3, comm, &send_request[2]);
        // receive from y-left proc
        MPI_Irecv(rbuff_y_left, nnode_x * nnode_z, MPI_DOUBLE, y_left, 2, comm, &recv_request[2]);
    }

    // exchange data with the right in y direction =================================
    if (y_right >= 0) {
        // collect data for sending
        count = 0;
        j = nnode_y;
        for (int i = 1; i <= nnode_x; i++) {
            for (int k = 1; k <= nnode_z; k++) {
                sbuff_y_right[count] = u[ID(i,j,k)];
                count = count + 1;
            }
        }
        // send to y-right proc
        MPI_Isend(sbuff_y_right, nnode_x * nnode_z, MPI_DOUBLE, y_right, 2, comm, &send_request[3]);
        // receive from y-right proc
        MPI_Irecv(rbuff_y_right, nnode_x * nnode_z, MPI_DOUBLE, y_right, 3, comm, &recv_request[3]);
    }

    // exchange data with the left in z direction=================================
    if (z_left >= 0) {
        // collect data for sending
        count = 0;
        k = 1;
        for (int i = 1; i <= nnode_x; i++) {
            for (int j = 1; j <= nnode_y; j++) {
                sbuff_z_left[count] = u[ID(i,j,k)];
                count = count + 1;
            }
        }
        // send to z-left proc
        MPI_Isend(sbuff_z_left, nnode_x * nnode_y, MPI_DOUBLE, z_left, 5, comm, &send_request[4]);
        // receive from z-left proc
        MPI_Irecv(rbuff_z_left, nnode_x * nnode_y, MPI_DOUBLE, z_left, 4, comm, &recv_request[4]);
    }

    // exchange data with the right in z direction =================================
    if (z_right >= 0) {
        // collect data for sending
        count = 0;
        k = nnode_z;
        for (int i = 1; i <= nnode_x; i++) {
            for (int j = 1; j <= nnode_y; j++) {
                sbuff_z_right[count] = u[ID(i,j,k)];
                count = count + 1;
            }
        }
        // send to z-right proc
        MPI_Isend(sbuff_z_right, nnode_x * nnode_y, MPI_DOUBLE, z_right, 4, comm, &send_request[5]);
        // receive from z-right proc
        MPI_Irecv(rbuff_z_right, nnode_x * nnode_y, MPI_DOUBLE, z_right, 5, comm, &recv_request[5]);
    }


    const double fac_1by6=1.0/6.0;
    // update displacement of interior nodes =================================
    for (int i = 2; i <= (nnode_x-1 ); i++) {
        for (int j = 2; j <= (nnode_y-1); j++) {
            for (int k = 2; k <= (nnode_z-1); k++) {
                // determine forcing
                x = (double)((i - 1) + ixmin)*h;
                y = (double)((j - 1) + iymin)*h;
                z = (double)((k - 1) + izmin)*h;
                b = sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
                // update u
                u_new[ID(i,j,k)] = fac_1by6 * (u[ID(i - 1,j,k)] + u[ID(i + 1,j,k)] +
                                            u[ID(i,j - 1,k)] + u[ID(i,j + 1,k)] +
                                           u[ID(i,j,k - 1)] + u[ID(i,j,k + 1)] - (h*h)*b);
                //printf("at center, x= %f, y= %f, z= %f\n",x, y, z);
            }
        }
    }



    // wait for sending and receive to be finished =================================
    if (x_left >= 0) {
        ierr = MPI_Wait(&recv_request[0],&status);
        ierr = MPI_Wait(&send_request[0],&status);
        // distribute received data (from left proc) to ghost nodes:
        count = 0;
        i = 0;          // data is put on the left side
        for (int j = 1; j <= nnode_y; j++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = rbuff_x_left[count];
                count = count + 1;
            }
        }
    }
    if (x_right >= 0) {
        ierr = MPI_Wait(&recv_request[1],&status);
        ierr = MPI_Wait(&send_request[1],&status);
        // distribute received data to ghost nodes:
        count = 0;
        i = nnode_x + 1; // data is put on the right side
        for (int j = 1; j <= nnode_y; j++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = rbuff_x_right[count];
                count = count + 1;
            }
        }
    }


    if (y_left >= 0) {
        ierr = MPI_Wait(&recv_request[2],&status);
        ierr = MPI_Wait(&send_request[2],&status);
        // distribute received data (from left proc) to ghost nodes:
        count = 0;
        j = 0;          // data is put on the left side
        for (int i = 1; i <= nnode_x; i++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = rbuff_y_left[count];
                count = count + 1;
            }
        }
    }
    if (y_right >= 0) {
        ierr = MPI_Wait(&recv_request[3],&status);
        ierr = MPI_Wait(&send_request[3],&status);
        // distribute received data (from left proc) to ghost nodes:
        count = 0;
        j = nnode_y + 1;          // data is put on the left side
        for (int i = 1; i <= nnode_x; i++) {
            for (int k = 1; k <= nnode_z; k++) {
                u[ID(i,j,k)] = rbuff_y_right[count];
                count = count + 1;
            }
        }
    }


    if (z_left >= 0) {
        ierr = MPI_Wait(&recv_request[4],&status);
        ierr = MPI_Wait(&send_request[4],&status);
        // distribute received data (from left proc) to ghost nodes:
        count = 0;
        k = 0;          // data is put on the left side
        for (int i = 1; i <= nnode_x; i++) {
            for (int j = 1; j <= nnode_y; j++) {
                u[ID(i,j,k)] = rbuff_z_left[count];
                count = count + 1;
            }
        }
    }
    if (z_right >= 0) {
        ierr = MPI_Wait(&recv_request[5],&status);
        ierr = MPI_Wait(&send_request[5],&status);
        // distribute received data (from left proc) to ghost nodes:
        count = 0;
        k = nnode_z + 1;          // data is put on the left side
        for (int i = 1; i <= nnode_x; i++) {
            for (int j = 1; j <= nnode_y; j++) {
                u[ID(i,j,k)] = rbuff_z_right[count];
                count = count + 1;
            }
        }
    }


    //CHECKING
    //if (rank == 2) {
    //    for (unsigned i=0; i < nnode_x*nnode_y*nnode_z; i++) {
    //        printf("rank= %d, u[%d]= %f\n",  rank, i, u[i]);
    //    }
    //}



    // now update boundary nodes =================================
    // x left side:

    i = 1;
    for (int j = 1; j <= nnode_y; j++) {
        for (int k = 1; k <= nnode_z; k++) {
            // determine forcing
            x = (double)((i - 1) + ixmin)*h;
            y = (double)((j - 1) + iymin)*h;
            z = (double)((k - 1) + izmin)*h;
            b = sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
            // update u
            u_new[ID(i,j,k)] = fac_1by6 * (u[ID(i - 1,j,k)] + u[ID(i + 1,j,k)] +
                                        u[ID(i,j - 1,k)] + u[ID(i,j + 1,k)] +
                                        u[ID(i,j,k - 1)] + u[ID(i,j,k + 1)] - (h*h)*b);
            //printf("at left %d , x= %f, y= %f, z= %f\n",rank,x, y, z);
        }
    }

    // x right side:
    i = nnode_x;
    for (int j = 1; j <= nnode_y; j++) {
        for (int k = 1; k <= nnode_z; k++) {
            // determine forcing
            x = (double)((i - 1) + ixmin)*h;
            y = (double)((j - 1) + iymin)*h;
            z = (double)((k - 1) + izmin)*h;
            b = sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
            // update u
            u_new[ID(i, j, k)] = fac_1by6 * (u[ID(i - 1,j,k)] + u[ID(i + 1,j,k)] +
                                          u[ID(i,j - 1,k)] + u[ID(i,j + 1,k)] +
                                          u[ID(i,j,k - 1)] + u[ID(i,j,k + 1)] - (h*h)*b);
            //printf("at center, x= %f, y= %f, z= %f\n",x, y, z);
        }
    }

    // y left side:
    j = 1;
    for (int i = 1; i <= nnode_x; i++) {
        for (int k = 1; k <= nnode_z; k++) {
            // determine forcing
            x = (double)((i - 1) + ixmin)*h;
            y = (double)((j - 1) + iymin)*h;
            z = (double)((k - 1) + izmin)*h;
            b = sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
            // update u
            u_new[ID(i, j, k)] = fac_1by6 * (u[ID(i - 1,j,k)] + u[ID(i + 1,j,k)] +
                                          u[ID(i,j - 1,k)] + u[ID(i,j + 1,k)] +
                                          u[ID(i,j,k - 1)] + u[ID(i,j,k + 1)] - (h*h)*b);
        }
    }

    // y right side:
    j = nnode_y;
    for (int i = 1; i <= nnode_x; i++) {
        for (int k = 1; k <= nnode_z; k++) {
            // determine forcing
            x = (double)((i - 1) + ixmin)*h;
            y = (double)((j - 1) + iymin)*h;
            z = (double)((k - 1) + izmin)*h;
            b = sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
            // update u
            u_new[ID(i, j, k)] = fac_1by6 * (u[ID(i - 1,j,k)] + u[ID(i + 1,j,k)] +
                                          u[ID(i,j - 1,k)] + u[ID(i,j + 1,k)] +
                                          u[ID(i,j,k - 1)] + u[ID(i,j,k + 1)] - (h*h)*b);
        }
    }

    // z left side:
    k = 1;
    for (int i = 1; i <= nnode_x; i++) {
        for (int j = 1; j <= nnode_y; j++) {
            // determine forcing
            x = (double)((i - 1) + ixmin)*h;
            y = (double)((j - 1) + iymin)*h;
            z = (double)((k - 1) + izmin)*h;
            b = sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
            // update u
            u_new[ID(i, j, k)] = fac_1by6 * (u[ID(i - 1,j,k)] + u[ID(i + 1,j,k)] +
                                          u[ID(i,j - 1,k)] + u[ID(i,j + 1,k)] +
                                          u[ID(i,j,k - 1)] + u[ID(i,j,k + 1)] - (h*h)*b);
        }
    }

    // z right side:
    k = nnode_z;
    for (int i = 1; i <= nnode_x; i++) {
        for (int j = 1; j <= nnode_y; j++) {
            // determine forcing
            x = (double)((i - 1) + ixmin)*h;
            y = (double)((j - 1) + iymin)*h;
            z = (double)((k - 1) + izmin)*h;
            b = sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
            // update u
            u_new[ID(i, j, k)] = fac_1by6 * (u[ID(i - 1,j,k)] + u[ID(i + 1,j,k)] +
                                          u[ID(i,j - 1,k)] + u[ID(i,j + 1,k)] +
                                          u[ID(i,j,k - 1)] + u[ID(i,j,k + 1)] - (h*h)*b);
        }
    }


        // apply Dirichlet boundary condition on u_new  =================================
        // (u_new now is not correct on the boundary because we also update the boundary points in the above calculation)
        // this update is not needed for the problem solving, but it is necessary for calculating the error |u_new - u| correctly

        if (x_left < 0) {
            i = 1;
            for (int j = 1; j <= nnode_y; j++) {
                for (int k = 1; k <= nnode_z; k++) {
                    u_new[ID(i,j,k)] = 0.0;
                }
            }
        }
        if (x_right < 0) {
            i = nnode_x;
            for (int j = 1; j <= nnode_y; j++) {
                for (int k = 1; k <= nnode_z; k++) {
                    u_new[ID(i,j,k)] = 0.0;
                }
            }
        }
        if (y_left < 0) {
            j = 1;
            for (int i = 1; i <= nnode_x; i++) {
                for (int k = 1; k <= nnode_z; k++) {
                    u_new[ID(i,j,k)] = 0.0;
                }
            }
        }
        if (y_right < 0) {
            j = nnode_y;
            for (int i = 1; i <= nnode_x; i++) {
                for (int k = 1; k <= nnode_z; k++) {
                    u_new[ID(i,j,k)] = 0.0;
                }
            }
        }
        if (z_left < 0) {
            k = 1;
            for (int i = 1; i <= nnode_x; i++) {
                for (int j = 1; j <= nnode_y; j++) {
                    u_new[ID(i,j,k)] = 0.0;
                }
            }
        }
        if (z_right < 0) {
            k = nnode_z;
            for (int i = 1; i <= nnode_x; i++) {
                for (int j = 1; j <= nnode_y; j++) {
                    u_new[ID(i,j,k)] = 0.0;
                }
            }
        }

    // Compute the relative error (using L2 norm)
    my_change = 0.0;
    for (int i = 1; i <= nnode_x; i++) {
        for (int j = 1; j <= nnode_y; j++) {
            for (int k = 1; k <= nnode_z; k++) {
                my_change += (u_new[ID(i,j,k)] - u[ID(i,j,k)]) * (u_new[ID(i,j,k)] - u[ID(i,j,k)]);
            }
        }
    }


    MPI_Allreduce (&my_change, &change, 1, MPI_DOUBLE, MPI_SUM, comm);
    change = sqrt(change);

    //if (rank == 0) {
    //    if(step%100==0)printf("step= %d, max of norm(u)= %f\n", step, change);
    //}

    // swap
    swap = u;
    u = u_new;
    u_new = swap;
    } while (epsilon < change && step < MAX_ITER);

    // Timing
    t2 = MPI_Wtime() - t1;
    MPI_Reduce((double*)&t2, (double*)&tmax, 1, MPI_DOUBLE, MPI_MAX,0,comm);
    if (rank == 0) {
        printf("Problem size N= %d; # processes= %d; # steps= %d, relative error= %.4f, time taken= %.6f\n",
                N, size, step, change, tmax);
    }

    free(u);
    free(u_new);

    free(sbuff_x_left);
    free(sbuff_x_right);
    free(rbuff_x_left);
    free(rbuff_x_right);

    free(sbuff_y_left);
    free(sbuff_y_right);
    free(rbuff_y_left);
    free(rbuff_y_right);

    free(sbuff_z_left);
    free(sbuff_z_right);
    free(rbuff_z_left);
    free(rbuff_z_right);

    MPI_Finalize();
    return 0;

}
