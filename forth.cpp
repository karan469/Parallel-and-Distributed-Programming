#include <bits/stdc++.h>
#include <omp.h>
#include<stdlib.h>
#include <time.h>
#include <chrono>

using namespace std;

double **make2dmatrix(long n);
void free2dmatrix(double ** M, long n);
void printmatrix(double **A, long n);

long matrix_size;
char algo;

void decomposeOpenMP(double **A, double **l, double **u, long n){
	cout<<"Hello";
}

void initializeVersion1(double **A, long n)
{
	srand48((unsigned int)time(NULL));
	long i, j;
	for (i=0;i<n;i++){
		for (j=0;j<n;j++){
            A[i][j] = drand48();
		}
	}
}

void initializeVersion2(double **A,long n){
	srand48((unsigned int)time(NULL));
	long i,j, k;
	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			 if(i==j){A[i][j]=1;}
            else if(j>i){A[i][j]=0;}
            else{A[i][j]=drand48();}
		}
	}
}

void initializeVersion3(double **A,long n){
	srand48((unsigned int)time(NULL));
	long i,j, k;
	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			 if(i==j){A[i][j] = drand48();}
            else if(j>i){A[i][j] = drand48();}
            else{A[i][j]=0;}
		}
	}
}


double **getMatrix(long size,int version)
{
	double **m=make2dmatrix(size); //allocated
	switch(version){
	case 1: //A
		initializeVersion1(m,size);
		break;
	case 2: //L
		initializeVersion2(m,size);
		break;
	case 3: //U
		initializeVersion3(m, size);
	default:
		printf("INVALID VERSION NUMBERee\n");
		exit(0);
	}
	return m;
}


int main(int argc, char** argv){
    srand48((unsigned int)time(NULL));
    matrix_size = strtol(argv[1], NULL, 10);
    long thread_count = strtol(argv[2], NULL, 10);
    if(thread_count<1){
		thread_count=5;
	}
    omp_set_num_threads(thread_count);
    // printf("%d\n", matrix_size);
    
    double **matrix=getMatrix(matrix_size,1);
	// double **l = getMatrix(matrix_size, 2);
	// double **u = getMatrix(matrix_size, 3);

    clock_t begin, end;
	double duration;
	begin = clock();

	// decomposeOpenMP(matrix,l, u, matrix_size);
	printmatrix(matrix, matrix_size);

	end = clock();
	duration = ((double)(end - begin)) / CLOCKS_PER_SEC;

    printf("\n**********************************\n\n");
	printf("Algo selected :%s\n","OpenMP");
	printf("Size of Matrix :%lu \n",matrix_size);
	printf("Number of Procs : %lu\n",thread_count);
	// printf("%s",check(matrix,matrix_size,version)==1? "DECOMPOSE SUCCESSFULL\n":"DECOMPOSE FAIL\n");
	printf("DECOMPOSE TIME TAKEN : %f seconds\n",duration);
	printf("\n**********************************\n\n");

    free2dmatrix(matrix,matrix_size);
	// free2dmatrix(l, matrix_size);
	// free2dmatrix(u, matrix_size);
    return 0;
}

double **make2dmatrix(long n)
{
	long i;
	double **m;
	m = (double**)malloc(n*sizeof(double*));
	for (i=0;i<n;i++)
		m[i] = (double*)malloc(n*sizeof(double));
	return m;
}

// only works for dynamic arrays:
void printmatrix(double **A, long n)
{
	printf("\n *************** MATRIX ****************\n\n");
	long i, j;
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
			printf("%f ",A[i][j]);
		printf("\n");
	}
}

void free2dmatrix(double ** M, long n)
{
	long i;
	if (!M) return;
	for(i=0;i<n;i++)
		free(M[i]);
	free(M);
}