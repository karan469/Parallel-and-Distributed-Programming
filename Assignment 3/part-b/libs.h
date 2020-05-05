#ifndef LIBS_H
#define LIBS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <limits.h>
#include <stdbool.h>
#include <assert.h>
#define comm MPI_COMM_WORLD
#define BUFSIZE INT_MAX

//Function for Matrix Multiplication
void Matrix_Multiply(float *A, float *B, float *C, int m, int n, int p){
	for (int i = 0; i < m; ++i)
	{
		for(int j = 0;j<p;j++){
			C[i*p+j] = 0;
			for(int k=0;k<n;k++){
				C[i*p+j] += A[i*n+k] * B[k*p + j];
			}
		}
	}
}

//Prints Matrix in 2 dimensional format
void printMatrix(float *A, int m, int n){
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			printf("%.5f ", A[j + i*n]);
		}
		printf("\n");
	}
}

//Adds two matrices
void addMatrices(float *A, float *B, float *C, int n){
	for(int i=0;i<n;i++){
		C[i] = A[i] + B[i];
	}
}

//Checks if serial and parallel answers are the same
bool isEqual(float *A, float *B, int n){
	for(int i=0;i<n;i++){
		if((int)(A[i]*100)!=(int)(B[i]*100)){
			return false;
		}
	}
	return true;
}

//Allocate Dynamic Memory for matrix
float* malloc_matrix(int m, int n){
	return (float*)malloc(m*n*sizeof(float));
}

#endif
