#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <random>

int main() {
const unsigned long  NUM=100000000;
int * A = (int *)malloc(sizeof(int)*NUM);

double t1= omp_get_wtime();

std::cout << "==== Outside #pragma omp parallel ====\n";
std::cout << "Number of threads can be used = " << omp_get_max_threads() << "\n";
std::cout << "Number of threads currently being used = " << omp_get_num_threads() << "\n";
std::cout << "\n";
std::cout << "Start paralleling...\n";

#pragma omp parallel
{
    std::cout << "Number of threads can be used = " << omp_get_max_threads() << "\n";
    std::cout << "Number of threads currently being used = " << omp_get_num_threads() << "\n";
    std::cout << "This is thread " << omp_get_thread_num() << "\n";

    int rank;
    const int max = omp_get_max_threads();
    rank = omp_get_thread_num(); // thread id

    // distribute the array  A to threads
    const unsigned long ib = (NUM*rank)/max; // beginning index for thread #rank
    const unsigned long ie =  (NUM*(rank+1))/max; // ending index for thread #rank

    for (unsigned long i=ib; i<ie; i++) {
        A[i] = rand() % 100; // generate random integers (0-99)
        //A[i] = 1;
    }
}
double t2=omp_get_wtime();

printf("Time taken: %.6fs\n", (t2-t1));
delete [] A;
return 0;
}
