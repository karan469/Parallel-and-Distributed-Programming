#ifndef LIBS_H
#define LIBS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <limits.h>
#define comm MPI_COMM_WORLD

void printMatrix(float *A, int m, int n){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			printf("%.1f ", A[j + i*n]);
		}
		printf("\n");
	}
}

float* malloc_matrix(int m, int n){
	return (float*)malloc(m*n*sizeof(float));
}

#endif
