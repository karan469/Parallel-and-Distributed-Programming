#include "stdc++.h"
#include <pthread.h>
#include <omp.h>
#include<stdlib.h>
#include <time.h>
#include <chrono>

// #define matrix_size 3
using namespace std;

int step_i=-1;
int NUM_OF_THREADS;
int matrix_size;
double **A,**l,**u,**origMatrix;
int *pi;
long n;


struct thread_data{
int loop_var;
int swap_var;
int core;
};

double **make2dmatrix();
void free2dmatrix(double ** M);
void printmatrix(double **A, long size);
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
		
		u[k][k] = colmax;
		
		for(long i=k+1;i<n;i++){
			l[i][k] = (double)(A[i][k]/colmax); //u(k,k) = max actually
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
	// int core = step_i++;
	// printf("%d\n",core);
	struct thread_data *t = (struct thread_data *) arg;

	// printf("%d %d\n",t->loop_var,t->swap_var);

	// printf("%d %d",core*(n/NUM_OF_THREADS),(core+1)*(n/NUM_OF_THREADS));

	// if(t->loop_var==1){
	// 	// cout << (t->core*n)/(NUM_OF_THREADS) << " " << ((t->core+1)*n)/(NUM_OF_THREADS) << endl;
	// }

	// cout <<t->loop_var << ":" << t->core << endl;

	for(long i=((t->core*n)/(NUM_OF_THREADS));i<(((t->core+1)*n)/(NUM_OF_THREADS));i++){
		double temp = A[t->loop_var][i];
		A[t->loop_var][i] = A[t->swap_var][i];
		A[t->swap_var][i] = temp;
	}

	for(long i=(t->core*t->loop_var)/(NUM_OF_THREADS);i<((t->core+1)*t->loop_var)/(NUM_OF_THREADS);i++){
		double gte = l[t->loop_var][i];
		l[t->loop_var][i] = l[t->swap_var][i];
		l[t->swap_var][i] = gte;
	}
	// cout << "A" << endl;
	// printmatrix(A,n);
	// cout << "L" << endl;
	// printmatrix(l,n);
	// cout << "U" << endl;
	// printmatrix(u,n);
	// cout<<"PI MATRIX: \n";
	// for(long i=0;i<n;i++){
	// 	cout<<pi[i]<<" ";
	// }

	pthread_exit(NULL);
		
}

void *multi2(void *arg){
	// int core = step_i++;
	// printf("%d\n",core);
	struct thread_data *t = (struct thread_data *) arg;

	for(long i=((t->loop_var+1) + (((n-t->loop_var-1)*t->core)/NUM_OF_THREADS));i<((t->loop_var+1) + (((n-t->loop_var-1)*(t->core+1))/NUM_OF_THREADS));i++){
		l[i][t->loop_var] = A[i][t->loop_var]/u[t->loop_var][t->loop_var]; //u(k,k) = max actually
		u[t->loop_var][i] = A[t->loop_var][i];
	}

	pthread_exit(NULL);
		
}


void *multi3(void *arg){
	// int core = step_i++;
	// printf("%d\n",core);
	struct thread_data *t = (struct thread_data *) arg;

	for(long i=((t->loop_var+1) + (((n-t->loop_var-1)*t->core)/NUM_OF_THREADS));i<((t->loop_var+1) + (((n-t->loop_var-1)*(t->core+1))/NUM_OF_THREADS));i++){
		// #pragma omp for
		for(long j=t->loop_var+1;j<n;j++){
			A[i][j] = A[i][j] - (l[i][t->loop_var])*(u[t->loop_var][j]);
		}
	}

	pthread_exit(NULL);
		
}	

void decomposePthread(/*double **A, double **l, double **u, int *pi, long n*/){
    pthread_t threads[NUM_OF_THREADS];
	struct thread_data td[NUM_OF_THREADS];

	
	cout<<"Hello\n";

	// printmatrix(A, n);
	// printmatrix(l, n);
	// printmatrix(u, n);


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

		// if(colmax==0){cerr<<"Singular Matrix\n";}
		
		//swap
		int temp = pi[k];
		pi[k] = pi[k_];
		pi[k_] = temp;

		// printmatrix(l,n);
		// printmatrix(u,n);

		// cout << "Reached " << endl;

		// printmatrix(A, n);

		int core = -1;
			
		for(int i=0;i<NUM_OF_THREADS;i++){
			td[i].loop_var = k;
			td[i].swap_var = k_;
			core++;
			td[i].core = core;
			// cout << td.core << endl;

			int rc = pthread_create(&threads[i], NULL, multi, (void *)&td[i]);

			if(rc){
				cerr << "Threads not created " << endl;
			}
		}

		//Joining and waiting for all threads to execute and complete
		for(int i=0;i<NUM_OF_THREADS;i++){
			pthread_join(threads[i],NULL);
		}
		// cerr << "A matrix " << endl;
		// printmatrix(A,n);
		// cerr << "L matrix " << endl;
		// printmatrix(l,n);
		// cerr << "U matrix " << endl;
		// printmatrix(u,n);

		u[k][k] = A[k][k];

		core = -1;
			
		for(int i=0;i<NUM_OF_THREADS;i++){
			td[i].loop_var = k;
			td[i].swap_var = k_;
			core++;
			td[i].core = core;
			// cout << td.core << endl;

			int rc = pthread_create(&threads[i], NULL, multi2, (void *)&td[i]);

			if(rc){
				cerr << "Threads not created " << endl;
			}
		}

		for(int i=0;i<NUM_OF_THREADS;i++){
			pthread_join(threads[i],NULL);
		}
		
		core = -1;
			
		for(int i=0;i<NUM_OF_THREADS;i++){
			td[i].loop_var = k;
			td[i].swap_var = k_;
			core++;
			td[i].core = core;
			// cout << td.core << endl;

			int rc = pthread_create(&threads[i], NULL, multi3, (void *)&td[i]);

			if(rc){
				cerr << "Threads not created " << endl;
			}
		}

		for(int i=0;i<NUM_OF_THREADS;i++){
			pthread_join(threads[i],NULL);
		}

		// for(long i=k+1;i<n;i++){
		// 	l[i][k] = A[i][k]/u[k][k]; //u(k,k) = max actually
		// 	u[k][i] = A[k][i];
		// }
		
		// for(long i=k+1;i<n;i++){
		// 	// #pragma omp for
		// 	for(long j=k+1;j<n;j++){
		// 		A[i][j] = A[i][j] - (l[i][k])*(u[k][j]);
		// 	}
		// }

	}

	// return td;
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

//initialize A matrix
void initializeVersion1(double **a, long n){
	srand48((unsigned int)time(NULL));
	long i, j;
	for (i=0;i<n;i++){
		for (j=0;j<n;j++){
            a[i][j] = abs(100*drand48());
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
void initializeVersion2(double **a,long n){
	// srand48((unsigned int)time(NULL));
	long i,j, k;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i==j){a[i][j]=1;}
            else{a[i][j]=0;}
            // else{A[i][j]=0;} //make 0
		}
	}
}

//initialize U matrix
void initializeVersion3(double **a,long n){
	// srand48((unsigned int)time(NULL));
	long i,j, k;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			// if(i<=j){A[i][j] = abs(100*drand48());}
            // else{A[i][j]=0;}
			a[i][j] = 0;
		}
	}
}

//initialize pi  matrix
void initializeVersion4(int *a, long n){
	long i;
	for(i=0;i<n;i++){
		a[i] = i;
	}
}

//get matrix A/L/U/pi
double **getMatrix(long size,int version)
{
	double **m=make2dmatrix(); //allocated
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

void checkAns(/*int *pi, double **A, double **l, double **u*/long size){
	
	double **P = make2dmatrix();
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

	double **mult1 = make2dmatrix();
	for(long i = 0; i < size; i++)
        for(long j = 0; j < size; j++)
            for(long k = 0; k < size; k++)
            {
                mult1[i][j] += P[i][k] * origMatrix[k][j];
            }
	
	// cout<<"\nP*A matrix ->\n";
	// printmatrix(mult1, size);

	double **mult2 = make2dmatrix();
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
	free2dmatrix(mult1);
	free2dmatrix(mult2);
	cout<<"Check sum: "<<gsum<<"\n";
}

int main(int argc, char** argv){
    srand48((unsigned int)time(NULL));
    n = strtol(argv[1], NULL, 10);
    int thread_count = strtol(argv[2], NULL, 10);
	NUM_OF_THREADS = thread_count;
	// printf("Thread count is %d",NUM_OF_THREADS);
    int rc;
	// struct thread_data td;

    ifstream fin;
    fin.open(argv[3]);

    int user_file = -1;
    if(fin){user_file = 0;}

    double **AA;
    AA = getMatrix(matrix_size, 1);

    for(int i=0;i<matrix_size;i++){
    	for(int j=0;j<matrix_size;j++){
    		fin>>AA[i][j];
    	}
    }

    fin.close();

    if(thread_count<1){
		thread_count=5;
	}


    // omp_set_num_threads(thread_count);
    // printf("%d\n", matrix_size);

    //initialization of matrices
	A=getMatrix(n,1);

	if(user_file==0){
		for(long p=0;p<matrix_size;p++){
			for(long w=0;w<matrix_size;w++){
				A[p][w] = AA[p][w];
			}
		}
	}

	// A = hardMatrix();
	// printmatrix(A,n);
	l = getMatrix(n, 2);
	// printmatrix(l,n);
	u = getMatrix(n, 3);
	// printmatrix(u,n);

	// double **checker = vec2Matrix();
	
	//Original A matrix stored
	origMatrix = getMatrix(n, 1);
	//recomended - dont use this extra loop for copying
	for(long p=0;p<n;p++){
		for(long w=0;w<n;w++){
			// origMatrix[p][w] = matrix[p][w];
			origMatrix[p][w] = A[p][w];

		}
	}
	
	// int *pi;
	pi = (int*)malloc(n*sizeof(int*));
	pi = (int*)malloc(n*sizeof(int));
	initializeVersion4(pi, n);

    // clock_t begin, end;
	double duration;
	double start = omp_get_wtime();

	//main function starts here
	// decomposePthread(/*A,l, u, pi, matrix_size*/);
	decomposePthread();
	// decomposeOpenMP(checker,l, u, pi, matrix_size);

	// // serialDecompose(matrix,l, u, pi, matrix_size);
	// printmatrix(origMatrix, n);
	// cerr<<"\nL Matrix:\n";
	// printmatrix(l, n);
	// cerr<<"\nU Matrix:\n";
	// printmatrix(u, n);

	double end = omp_get_wtime(); 
	duration = ((double)(end - start));

	//Printing results
    printf("\n**********************************\n\n");
	printf("Algo selected :%s\n","Pthread");
	printf("Size of Matrix :%lu \n",n);
	printf("Number of Procs : %lu\n",NUM_OF_THREADS);
	// printf("%s",check(matrix,matrix_size,version)==1? "DECOMPOSE SUCCESSFULL\n":"DECOMPOSE FAIL\n");
	printf("DECOMPOSE TIME TAKEN : %f seconds\n",duration);
	printf("\n**********************************\n\n");

	//Checking the value of residual matrix i.e. L2,1 norm of (PA-LU)
	checkAns(n);

	if(user_file==0){
		ofstream fout;
		fout.open("./dump/P_" + to_string(matrix_size) + "_" + to_string(thread_count) + ".txt");

		for(int i = 0;i<matrix_size;i++){
			fout<<pi[i]<<" ";
		}

		fout.close();

		fout.open("./dump/L_" + to_string(matrix_size) + "_" + to_string(thread_count) + ".txt");
		for(int i=0;i<matrix_size;i++){
			for(int j=0;j<matrix_size;j++){
				fout<<l[i][j]<<" ";
			}
			fout<<"\n";
		}

		fout.close();

		fout.open("./dump/U_" + to_string(matrix_size) + "_" + to_string(thread_count) + ".txt");
		for(int i=0;i<matrix_size;i++){
			for(int j=0;j<matrix_size;j++){
				fout<<u[i][j]<<" ";
			}
			fout<<"\n";
		}

		fout.close();
	}

	//Freeing 2dm atrices
    free2dmatrix(A);
	free2dmatrix(origMatrix);
    // free2dmatrix(checker,matrix_size);
	free2dmatrix(l);
	free2dmatrix(u);
    return 0;
}

//Allocate mem space for a 2-D matrix
double **make2dmatrix()
{
	long i;
	double **m;
	m = (double**)malloc(n*sizeof(double*));
	for (i=0;i<n;i++)
		m[i] = (double*)malloc(n*sizeof(double));
	return m;
}

// only works for dynamic arrays:
void printmatrix(double **M, long size)
{
	printf("\n *************** MATRIX ****************\n\n");
	long i, j;
	for (i=0;i<size;i++)
	{
		for (j=0;j<size;j++)
			printf("%f ",M[i][j]);
		printf("\n");
	}
}

void free2dmatrix(double ** M)
{
	long i;
	if (!M) return;
	for(i=0;i<n;i++)
		free(M[i]);
	free(M);
}
