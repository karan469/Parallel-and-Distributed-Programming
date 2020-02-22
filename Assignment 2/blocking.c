#include "libs.h"
#define BUFSIZE INT_MAX

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
	// A: (N x M) B: (M x N)
	const int M = 4;
	int N =  4;

	int rank, num_processes;
	MPI_Init(NULL, NULL);

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &num_processes);

	if (num_processes < 2) {
	    fprintf(stderr, "Must use atleast two processes for this example\n");
	    MPI_Abort(comm, 1);
	}

	// float buf[BUFSIZE];
	// MPI_Buffer_attach(buf, BUFSIZE);

	if(rank == 0)
	{
		int a_mssg_id = 0;
		int b_mssg_id = 0;

		A = malloc_matrix(N, M);
		B = malloc_matrix(M, N);
		C = malloc_matrix(N, N);

		srand(time(NULL));
		for(int a = 0;a<(N*M);a++){
			A[a] = rand()/(float)RAND_MAX * 10;
		}

		for(int b = 0;b<(M*N);b++){
			B[b] = rand()/(float)RAND_MAX * 10;
		}

		for(int c = 0;c<(M*M);c++){
			C[c] = (float)0;
		}

		printf("================== MATRIX A ===============\n");
		printMatrix(A, N, M);
		printf("================== MATRIX B ===============\n");
		printMatrix(B, M, N);
		// printMatrix(C, M, M);
		
		for(int f=0;f<num_processes-1;f++){
			A_block = malloc_matrix(N, M/(num_processes-1));
			B_block = malloc_matrix(M/(num_processes-1), N);

			for(int i=0;i<N*M/2;i++){
				B_block[i] = B[i + (N*M*f)/2];
			}

			int counter = 0;
			for(int i=0;i<N;i++){
				for(int j = i*M + f*(M/2); j<i*M + (int)(M/2) + f*(M/2);j++){
					A_block[counter++] = A[j];
				}
			}

			MPI_Rsend(A_block, (int)(N*M/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*13, comm);
			MPI_Rsend(B_block, (int)(M*N/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*97, comm);
		}
	}
	else if (rank!=0)
	{
		printf("process %d started..\n", rank);
		A_block = malloc_matrix(N/(num_processes-1), M);
		B_block = malloc_matrix(N/(num_processes-1), M);
		C_block = malloc_matrix(N/(num_processes-1), M);

		MPI_Status status;
		MPI_Recv(A_block, INT_MAX, MPI_FLOAT, 0, (rank)*13, comm, &status);
		MPI_Recv(B_block, INT_MAX, MPI_FLOAT, 0, (rank)*97, comm, &status);
		printf("====================== Process %d ===================\n", rank);
		printMatrix(A_block, N, M/(num_processes-1));
		printMatrix(B_block, M/(num_processes-1), N);
		printf("\n");

	}


	MPI_Finalize();
	return 0;
}

// 1:
// A -> 0
// B -> 1
// 2:
// A -> 2
// B -> 3
