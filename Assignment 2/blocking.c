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

		printMatrix(A, N, M);
		// printMatrix(B, M, N);
		// printMatrix(C, M, M);
		
		for(int f=0;f<num_processes-1;f++){
			// printf("YOLO\n");
			MPI_Rsend(A, (int)(N*M/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*13, comm);
			MPI_Rsend(B, (int)(M*N/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*97, comm);
			
			// int m_line = 0;
			// for(int k=0;k<N*M/(num_processes-1);k++){
			// 	MPI_Rsend(A[], )
			// }
			for(int k=0;k<N;k++){
				MPI_Rsend(A[f*(M/(num_processes-1)) + k*M], (int)(M/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*13, comm);
			}
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
		printMatrix(A_block, N/(num_processes-1), M);
		printMatrix(B_block, M, N/(num_processes-1));
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