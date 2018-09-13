#include <iostream>
#include <vector>
#include <functional>
#include <omp.h>
//#include <cmath>
#include <time.h>
#include "SpMv.h"
#include <fstream>
#include <string>
using namespace std;

int main() {

    unsigned int n, nnz;
    const unsigned int tol=0.0001;

    // Open the text file containing data of sparse matrix
    ifstream myfile;
    myfile.open("1138_bus.txt");
    if (!myfile ) {
        cout << "Unable to open file" << std::endl;
        exit(1);
    }
    // read the size of matrix and number of non-zero elements
    myfile >> n >> nnz;

    // allocate matrix and vectors
    double *val = (double *)malloc(sizeof(double)*nnz);
    unsigned int *col = (unsigned int *)malloc(sizeof(unsigned int)*nnz);
    unsigned int *rpt = (unsigned int *)malloc(sizeof(unsigned int)*(n+1));

    // read the row-index range
    for (int i=0; i < n+1; i++)
        myfile >> rpt[i];
    // read the column-index and value of non-zero elements
    for (int i=0; i < nnz; i++)
        myfile >> col[i] >> val[i];

    // convert row and column index to C/C++ style (i.e. starting from 0)
    for (int i=0; i < n+1; i++)
        rpt[i] = rpt[i] - 1;
    for (int i=0; i < nnz; i++)
        col[i] = col[i] - 1;


    // Generate random vector X
    double *X = (double *)malloc(sizeof(double)*n);
    for (unsigned int i=0; i<n; i++) {
        X[i] = (double)std::rand()/(double)(RAND_MAX/5.0);
    }
    // Result vector Y
    double *Y = (double *)malloc(sizeof(double)*n);

    // Matrix-vector multiplication
    printf("np = %d\n",omp_get_max_threads());
    double t1= omp_get_wtime();
    SpMv(val,col,rpt, X, Y, n);
    double t2=omp_get_wtime();
    printf("Time taken: %.6fs\n", (t2-t1));

    //==================== Testing ====================
    for (int i=0; i < n+1; i++)
        rpt[i] = i;
    for (int i=0; i < nnz; i++) {
        col[i] = i;
        val[i] = 1;
    }
    SpMv(val,col,rpt, X, Y, n);
    // difference between X and Y
    for (int i=0; i<n; i++)
        if (Y[i]-X[i] > tol) {
            cout << "The test does not pass.";
            exit(1);
        }
    printf("The test is passed!\n");
    //==================== End of Testing ====================

    delete [] val;
    delete [] col;
    delete [] rpt;
    delete [] X;
    delete [] Y;

    return 0;

}
