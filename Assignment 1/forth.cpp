#include <bits/stdc++.h>
#include <omp.h>
// #include<stdlib.h>
#include <time.h>
#include <chrono>
#include <fstream>

using namespace std;

double **make2dmatrix(long n);
void free2dmatrix(double ** M, long n);
void printmatrix(double **A, long n);

long matrix_size;
char algo;

void serialDecompose(double **A, double **l, double **u, int *pi, long n){
	// cout<<"Hello\n";
	
	for(int k=0;k<n;k++){
		long rows,colmax;
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

//Openmp parallel code for decomposing. Input - A(n, n) | Output - pi(n), L(n, n), U(n, n)
void decomposeOpenMP(double **A, double **l, double **u, int *pi, long n){
	double colmax;
	int k_;
	for(long k=0;k<n;k++){
		colmax = 0;

		for(int i=k;i<n;i++){
			if(colmax<abs(A[i][k])){colmax = A[i][k];k_ = i;}
		}
		if(colmax==0){cerr<<"Singular Matrix\n";}
		
		//swap
		// int temp = pi[k];
		// pi[k] = pi[k_];
		// pi[k_] = temp;
		swap(pi[k], pi[k_]);
		#pragma omp parallel default(none) shared(A, l, u, n, colmax, k, k_)
		{
			#pragma omp for nowait
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
			// #pragma omp barrier
			#pragma omp for
			for(long i=k;i<n;i++){
				// l[i][k] = A[i][k]/u[k][k];
				l[i][k] = (double)(A[i][k]/colmax);
				u[k][i] = A[k][i];
			}
			
			#pragma omp for collapse(2)
			for(long i=k+1;i<n;i++){
				for(long j=k+1;j<n;j++){
					A[i][j] -= l[i][k]*u[k][j];
				}
			}
		}
	}
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
		}
	}

	// cout<<"\nL*U Matrix ->\n";
	// printmatrix(mult2, size);

	for(long i = 0; i < size; i++)
    {    for(long j = 0; j < size; j++){
			mult1[i][j] -= mult2[i][j];
		}
	}
	double long gsum = 0.000000;
	for(long i=0;i<size;i++){
		double long sum = 0.000000;
		for(long j=0;j<size;j++){
			sum += mult1[i][j]*mult1[i][j];
		}
		gsum += sqrt(sum);
	}
	free2dmatrix(mult1,size);
	free2dmatrix(mult2,size);
	cout<<"Check sum: "<<gsum<<"\n";
	if(gsum < 4.5e-10){cout<<"DECOMPOSITION CORRECTLY DONE."<<endl;}
}

int main(int argc, char** argv){

    srand48((unsigned int)time(NULL));
    matrix_size = strtol(argv[1], NULL, 10);
    long thread_count = strtol(argv[2], NULL, 10);


    fstream fin;
    fin.open(argv[3]);

    int user_file = -1;
    if(fin){user_file = 0;}

    double **A;
    A = getMatrix(matrix_size, 1);

    for(int i=0;i<matrix_size;i++){
    	for(int j=0;j<matrix_size;j++){
    		fin>>A[i][j];
    	}
    }

    // cout<<A[0][99]<<endl;
    
    if(thread_count<1){
		thread_count=5;
	}
    omp_set_num_threads(thread_count);
    
    //initialization of matrices
	double **matrix=getMatrix(matrix_size,1);
	
	if(user_file==0){
		for(long p=0;p<matrix_size;p++){
			for(long w=0;w<matrix_size;w++){
				matrix[p][w] = A[p][w];
			}
		}
	}

	double **l = getMatrix(matrix_size, 2);
	double **u = getMatrix(matrix_size, 3);
	
	//Original A matrix stored
	double **origMatrix = getMatrix(matrix_size, 1);
	//recomended - dont use this extra loop for copying
	for(long p=0;p<matrix_size;p++){
		for(long w=0;w<matrix_size;w++){
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
	decomposeOpenMP(matrix,l, u, pi, matrix_size);
	// serialDecompose(matrix,l, u, pi, matrix_size);
	// printmatrix(matrix, matrix_size);
	// cerr<<"\nL Matrix:\n";
	// printmatrix(l, matrix_size);
	// cerr<<"\nU Matrix:\n";
	// printmatrix(u, matrix_size);

	double end = omp_get_wtime(); 
	duration = ((double)(end - start));

	//Printing results
    printf("\n**********************************\n\n");
	printf("Algo selected :%s\n","OpenMP");
	printf("Size of Matrix :%lu \n",matrix_size);
	printf("Number of Procs : %lu\n",thread_count);
	printf("DECOMPOSE TIME TAKEN : %f seconds\n",duration);
	printf("\n**********************************\n\n");

	//Checking the value of residual matrix i.e. L2,1 norm of (PA-LU)
	checkAns(pi, origMatrix, l, u, matrix_size);

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
    free2dmatrix(matrix,matrix_size);
	free2dmatrix(l, matrix_size);
	free2dmatrix(u, matrix_size);
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
