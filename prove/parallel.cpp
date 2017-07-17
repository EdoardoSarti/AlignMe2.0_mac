#include <omp.h>
#include <stdio.h>

int main() {
    int t, nt;
    omp_set_num_threads(10);
#pragma omp parallel
	printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
}
