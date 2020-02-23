#include "libs.h"

float *A, *B, *C;
float *A_block, *B_block, *C_block;

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

int main(int argc, char const *argv[])
{
	const int M = 32;
	
	int N =  atoi(argv[1]);
	

	int rank, num_processes;
	MPI_Init(NULL, NULL);

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &num_processes);

	if (num_processes < 2) {
	    fprintf(stderr, "Must use atleast two processes for this example\n");
	    MPI_Abort(comm, 1);
	}

	if(rank == 0)
	{
		int a_mssg_id = 0;
		int b_mssg_id = 0;

		A = malloc_matrix(N, M);
		B = malloc_matrix(M, N);
		C = malloc_matrix(N, N);

		srand(time(NULL));
		for(int a = 0;a<(N*M);a++){
			A[a] = rand()/(float)RAND_MAX;
		}

		for(int b = 0;b<(M*N);b++){
			B[b] = rand()/(float)RAND_MAX;
		}

		for(int c = 0;c<(N*N);c++){
			C[c] = (float)0;
		}

		// printf("================== MATRIX A ===============\n");
		// printMatrix(A, N, M);
		// printf("================== MATRIX B ===============\n");
		// printMatrix(B, M, N);
		// printMatrix(C, N, N);
		
		double start = MPI_Wtime();

		for(int f=0;f<num_processes-1;f++){
			A_block = malloc_matrix(N, M/(num_processes-1));
			B_block = malloc_matrix(M/(num_processes-1), N);

			for(int i=0;i<N*M/(num_processes-1);i++){
				B_block[i] = B[i + (N*M*f)/(num_processes-1)];
			}

			int counter = 0;
			for(int i=0;i<N;i++){
				for(int j = i*M + f*(M/(num_processes-1)); j<i*M + (int)(M/(num_processes-1)) + f*(M/(num_processes-1));j++){
					A_block[counter++] = A[j];
				}
			}

			MPI_Rsend(A_block, (int)(N*M/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*13, comm);
			MPI_Rsend(B_block, (int)(M*N/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*97, comm);
		}

		C_block = malloc_matrix(N, N);
		for(int c = 0;c<(N*N);c++){
			C_block[c] = (float)0;
		}

		MPI_Status status;
		for(int f=0;f<num_processes-1;f++){
			MPI_Recv(C_block, N*N, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
			addMatrices(C, C_block, C, N*N);
		}

		double end = MPI_Wtime();

        double duration = (float)end-start;

        printf("[BLOCKING P2P] Time it took to run the process is %0.4fs\n",duration);

		// float *D;
		// D = malloc_matrix(N, N);
		// for(int d = 0;d<(N*N);d++){
		// 	D[d] = (float)0;
		// }
		
		// Matrix_Multiply(A, B, D, N, M, N);
		// printf("%d\n", isEqual(C, D, N));
	}
	else if (rank>0)
	{
		A_block = malloc_matrix(N, M/(num_processes-1));
		B_block = malloc_matrix(M/(num_processes-1), N);
		C_block = malloc_matrix(N,N);

		MPI_Status status;
		MPI_Recv(A_block, INT_MAX, MPI_FLOAT, 0, (rank)*13, comm, &status);
		MPI_Recv(B_block, INT_MAX, MPI_FLOAT, 0, (rank)*97, comm, &status);

		for(int c = 0;c<(N*N);c++){
			C_block[c] = (float)0;
		}
		Matrix_Multiply(A_block, B_block, C_block, N, M/(num_processes-1), N);
		MPI_Send(C_block, N*N, MPI_FLOAT, 0, (rank)*113, comm);
	}


	MPI_Finalize();
	return 0;
}

