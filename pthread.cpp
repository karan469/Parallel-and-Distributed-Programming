#include "stdc++.h"
#include <pthread.h>
#include <omp.h>
#include<stdlib.h>
#include <time.h>
#include <chrono>

// #define matrix_size 3
using namespace std;

int step_i=0;
int NUM_OF_THREADS;
int matrix_size;

struct thread_data{
	double **A;
	double **l;
	double **u;
	int *pi;
	long n;
	long loop_var;
	long swap_var;
};

double **make2dmatrix(long n);
void free2dmatrix(double ** M, long n);
void printmatrix(double **A, long n);
void* multi(void* arg);

char algo;

void serialDecompose(double **A, double **l, double **u, int *pi, long n){
	cout<<"Hello\n";
	
	for(int k=0;k<n;k++){
		long rows;
		double colmax;
		int temp;
		int k_;
		colmax = 0;
		for(long i=k;i<n;i++){
			if(colmax<abs(A[i][k])){colmax = A[i][k];k_ = i;}
		}
		if(colmax==0){cerr<<"Singular Matrix\n";}
		
		//swap
		temp = pi[k];
		pi[k] = pi[k_];
		pi[k_] = temp;

		// cout<<"PI MATRIX "<<k<<"rd iteration: \n";
		// for(int i=0;i<n;i++){
		// 	cout<<pi[i]<<" ";
		// }
		// cout<<"\n";

		for(long i=0;i<n;i++){
			double temp = A[k][i];
			A[k][i] = A[k_][i];
			A[k_][i] = temp;
		}

		
		for(long i=0;i<k;i++){ //dwijesh
			// l[k][i] += l[k_][i];
			// l[k_][i] = l[k][i]-l[k_][i];
			// l[k][i] -= l[k_][i];
			double gte = l[k][i];
			l[k][i] = l[k_][i];
			l[k_][i] = gte;
		}
		
		u[k][k] = A[k][k];
		
		for(long i=k+1;i<n;i++){
			l[i][k] = A[i][k]/u[k][k]; //u(k,k) = max actually
			u[k][i] = A[k][i];
		}
		
		for(long i=k+1;i<n;i++){
			// #pragma omp for
			for(long j=k+1;j<n;j++){
				A[i][j] = A[i][j] - l[i][k]*u[k][j];
			}
		}	
	}
}

void *multi(void *arg){
	int core = step_i++;
	struct thread_data *t = (struct thread_data *) arg;

	for(long i=core*(matrix_size/NUM_OF_THREADS);i<(core+1)*(matrix_size/NUM_OF_THREADS);i++){
		double temp = t->A[t->loop_var][i];
		t->A[t->loop_var][i] = t->A[t->swap_var][i];
		t->A[t->swap_var][i] = temp;
	}

	for(long i=core*(t->loop_var/NUM_OF_THREADS);i<(core+1)*(t->loop_var/NUM_OF_THREADS);i++){
		double gte = t->l[t->loop_var][i];
		t->l[t->loop_var][i] = t->l[t->swap_var][i];
		t->l[t->swap_var][i] = gte;
	}

	pthread_exit(NULL);
		
}

struct thread_data decomposePthread(double **A, double **l, double **u, int *pi, long n){
    pthread_t threads[NUM_OF_THREADS];
	struct thread_data td;
	
	cout<<"Hello\n";
	// int k;
	for(long k=0;k<n;k++){
		long rows,mymin,mymax,colmax;
		// int pid=0;
		// int nprocs;
		// int temp;
		int k_;
		colmax = 0;
		for(int i=k;i<n;i++){
			if(colmax<abs(A[i][k])){colmax = A[i][k];k_ = i;}
		}
		if(colmax==0){cerr<<"Singular Matrix\n";}
		
		//swap
		int temp = pi[k];
		pi[k] = pi[k_];
		pi[k_] = temp;

		td.A = A;
		td.l = l;
		td.u = u;
		td.pi = pi;
		td.n = n;	
		td.loop_var = k;
		td.swap_var = k_;

		for(int i=0;i<NUM_OF_THREADS;i++){
			int rc = pthread_create(&threads[i], NULL, multi, (void *)&td);

			if(rc){
				cerr << "Threads not created " << endl;
			}
		}

		//Joining and waiting for all threads to execute and complete
		for(int i=0;i<NUM_OF_THREADS;i++){
			pthread_join(threads[i],NULL);
		}

		td.u[k][k] = td.A[k][k];
		
		for(long i=k+1;i<n;i++){
			td.l[i][k] = td.A[i][k]/td.u[k][k]; //u(k,k) = max actually
			td.u[k][i] = td.A[k][i];
		}
		
		for(long i=k+1;i<n;i++){
			// #pragma omp for
			for(long j=k+1;j<n;j++){
				td.A[i][j] = td.A[i][j] - (td.l[i][k])*(td.u[k][j]);
			}
		}

	}

	return td;
}

//Openmp parallel code for decomposing. Input - A(n, n) | Output - pi(n), L(n, n), U(n, n)
void decomposeOpenMP(double **A, double **l, double **u, int *pi, long n){
	cout<<"Hello\n";
	// int k;
	for(long k=0;k<n;k++){
		long rows,mymin,mymax,colmax;
		// int pid=0;
		// int nprocs;
		// int temp;
		int k_;
		colmax = 0;
		for(int i=k;i<n;i++){
			if(colmax<abs(A[i][k])){colmax = A[i][k];k_ = i;}
		}
		if(colmax==0){cerr<<"Singular Matrix\n";}
		
		//swap
		int temp = pi[k];
		pi[k] = pi[k_];
		pi[k_] = temp;

		#pragma omp parallel default(none) shared(A, l, u, n, colmax, k, k_)
		{
			// nprocs = omp_get_num_threads();
			// pid = omp_get_thread_num();
			// #pragma omp sections
			// {
			// 	#pragma omp section
			// 	{
			// 		#pragma omp for nowait //no wait can be implemented
			// 		for(long i=0;i<n;i++){
			// 			#pragma omp critical
			// 			{
			// 				double temp = A[k][i];
			// 				A[k][i] = A[k_][i];
			// 				A[k_][i] = temp;
			// 			}
			// 		}		
			// 	}

			// 	#pragma omp section
			// 	{
			// 		#pragma omp for
			// 		for(long i=0;i<k;i++){
			// 			double gte = l[k][i];
			// 			l[k][i] = l[k_][i];
			// 			l[k_][i] = gte;
			// 		}
			// 	}
			// }
			#pragma omp for nowait //no wait can be implemented
			for(long i=0;i<n;i++){
				#pragma omp critical
				{
					double temp = A[k][i];
					A[k][i] = A[k_][i];
					A[k_][i] = temp;
				}
			}

			#pragma omp for
			for(long i=0;i<k;i++){
				double gte = l[k][i];
				l[k][i] = l[k_][i];
				l[k_][i] = gte;
			}


			u[k][k] = colmax;
			
			#pragma omp for
			for(long i=k+1;i<n;i++){
				// l[i][k] = A[i][k]/u[k][k];
				l[i][k] = A[i][k]/colmax;
				u[k][i] = A[k][i];
			}
			
			#pragma omp for collapse(2)
			for(long i=k+1;i<n;i++){
				// #pragma omp for
				for(long j=k+1;j<n;j++){
					A[i][j] -= l[i][k]*u[k][j];
				}
			}
		}
	}
}

double **vec2Matrix(){
double **m;

	m = (double**)malloc(3*sizeof(double*));
	for (int i=0;i<3;i++)
		m[i] = (double*)malloc(3*sizeof(double));

	m[0][0]=10.000; 
	m[0][1]=-7.000; 
	m[0][2]=0.000; 
	m[1][0]=-3.000; 
	m[1][1]=2.000; 
	m[1][2]=6.000; 
	m[2][0]=5.000; 
	m[2][1]=-1.000; 
	m[2][2]=5.000; 

	return m;

}



//initialize A matrix
void initializeVersion1(double **A, long n){
	srand48((unsigned int)time(NULL));
	long i, j;
	for (i=0;i<n;i++){
		for (j=0;j<n;j++){
            A[i][j] = abs(100*drand48());
		}
	}
	// A[0][0] = 1;
	// A[0][1] = 2;
	// A[0][2] = 3;
	// A[1][0] = 3;
	// A[1][1] = 1;
	// A[1][2] = 2;
	// A[2][0] = 2;
	// A[2][1] = 3;
	// A[2][2] = 1;
}

//initialize L matrix
void initializeVersion2(double **A,long n){
	srand48((unsigned int)time(NULL));
	long i,j, k;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i==j){A[i][j]=1;}
            else{A[i][j]=0;}
            // else{A[i][j]=0;} //make 0 - dwijesh
		}
	}
}

//initialize U matrix
void initializeVersion3(double **A,long n){
	srand48((unsigned int)time(NULL));
	long i,j, k;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			// if(i<=j){A[i][j] = abs(100*drand48());}
            // else{A[i][j]=0;}
			A[i][j] = 0;
		}
	}
}

//initialize pi  matrix
void initializeVersion4(int *A, long n){
	long i;
	for(i=0;i<n;i++){
		A[i] = i;
	}
}

//get matrix A/L/U/pi
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
		break;
	default:
		printf("INVALID VERSION NUMBERee\n");
		exit(0);
	}
	return m;
}

void checkAns(int *pi, double **A, double **l, double **u, long size){
	double **P = make2dmatrix(size);
	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){
			P[i][j] = 0;
		}
		P[i][pi[i]] = 1;
	}
	// cout<<"PI MATRIX: \n";
	// for(int i=0;i<size;i++){
	// 	cout<<pi[i]<<" ";
	// }
	// cout<<"\n";

	// cerr<<"\nP matrix:\n";
	// printmatrix(P, size);

	double **mult1 = make2dmatrix(size);
	for(long i = 0; i < size; i++)
        for(long j = 0; j < size; j++)
            for(long k = 0; k < size; k++)
            {
                mult1[i][j] += P[i][k] * A[k][j];
            }
	
	// cout<<"\nP*A matrix ->\n";
	// printmatrix(mult1, size);

	double **mult2 = make2dmatrix(size);
	for(long i = 0; i < size; i++)
    {    for(long j = 0; j < size; j++)
        {    for(long k = 0; k < size; k++)
            {
                mult2[i][j] += l[i][k] * u[k][j];
            }
			// cout<<mult2[i][j]<<" ";
		}
		// cout<<"\n";
	}

	// cout<<"\nL*U Matrix ->\n";
	// printmatrix(mult2, size);

	for(long i = 0; i < size; i++)
    {    for(long j = 0; j < size; j++){
			mult1[i][j] -= mult2[i][j];
			// cout<<mult1[i][j]<<" ";
		}
		// cout<<"\n";
	}
	long  gsum = 0;
	for(long i=0;i<size;i++){
		long  sum = 0;
		for(long j=0;j<size;j++){
			sum += mult1[i][j]*mult1[i][j];
		}
		gsum += sqrt(sum);
		// cout << "Gsum is " << gsum << "\n";
	}
	free2dmatrix(mult1,size);
	free2dmatrix(mult2,size);
	cout<<"Check sum: "<<gsum<<"\n";
}

int main(int argc, char** argv){
    srand48((unsigned int)time(NULL));
    matrix_size = strtol(argv[1], NULL, 10);
    long thread_count = strtol(argv[2], NULL, 10);
	NUM_OF_THREADS = thread_count;
    int rc;
	struct thread_data td;

    if(thread_count<1){
		thread_count=5;
	}


    // omp_set_num_threads(thread_count);
    // printf("%d\n", matrix_size);

    //initialization of matrices
	double **matrix=getMatrix(matrix_size,1);
	double **l = getMatrix(matrix_size, 2);
	double **u = getMatrix(matrix_size, 3);
	// double **checker = vec2Matrix();
	
	//Original A matrix stored
	double **origMatrix = getMatrix(matrix_size, 1);
	//recomended - dont use this extra loop for copying
	for(long p=0;p<matrix_size;p++){
		for(long w=0;w<matrix_size;w++){
			// origMatrix[p][w] = matrix[p][w];
			origMatrix[p][w] = matrix[p][w];

		}
	}
	
	int *pi;
	pi = (int*)malloc(matrix_size*sizeof(int*));
	pi = (int*)malloc(matrix_size*sizeof(int));
	initializeVersion4(pi, matrix_size);

    // clock_t begin, end;
	double duration;
	double start = omp_get_wtime();

	//main function starts here
	td = decomposePthread(matrix,l, u, pi, matrix_size);
	// decomposeOpenMP(checker,l, u, pi, matrix_size);

	// serialDecompose(matrix,l, u, pi, matrix_size);
	// printmatrix(origMatrix, matrix_size);
	// cerr<<"\nL Matrix:\n";
	// printmatrix(td.l, matrix_size);
	// cerr<<"\nU Matrix:\n";
	// printmatrix(td.u, matrix_size);

	double end = omp_get_wtime(); 
	duration = ((double)(end - start));

	//Printing results
    printf("\n**********************************\n\n");
	printf("Algo selected :%s\n","OpenMP");
	printf("Size of Matrix :%lu \n",matrix_size);
	printf("Number of Procs : %lu\n",thread_count);
	// printf("%s",check(matrix,matrix_size,version)==1? "DECOMPOSE SUCCESSFULL\n":"DECOMPOSE FAIL\n");
	printf("DECOMPOSE TIME TAKEN : %f seconds\n",duration);
	printf("\n**********************************\n\n");

	//Checking the value of residual matrix i.e. L2,1 norm of (PA-LU)
	checkAns(td.pi, origMatrix, td.l, td.u, matrix_size);

	//Freeing 2dm atrices
    free2dmatrix(td.A,matrix_size);
    // free2dmatrix(checker,matrix_size);
	free2dmatrix(td.l, matrix_size);
	free2dmatrix(td.u, matrix_size);
    return 0;
}

//Allocate mem space for a 2-D matrix
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
