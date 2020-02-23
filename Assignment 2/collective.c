#include "libs.h"

float *A, *B, *C;
float *A_block, *B_block, *C_block;

void compute_ans(float *big, float *ans, int world_size, int N){
	for(int i=0;i<N*N;i++){
		ans[i] = (float)0;
	}

	for(int i=0;i<N*N;i++){
		float sum = 0;
		for(int j=i;j<world_size*N*N;j+=N*N){
			sum += big[j];
		}
		ans[i] = sum;
	}
}

int main(int argc, char const *argv[])
{
	const int M = 4;
	
	int N =  atoi(argv[1]);

	int rank, num_processes;
	MPI_Init(NULL, NULL);

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &num_processes);

	if (num_processes < 2) {
	    fprintf(stderr, "Must use atleast two processes for this example\n");
	    MPI_Abort(comm, 1);
	}


	B = malloc_matrix(M, N);
	A_block = malloc_matrix(N/(num_processes-1), M);
	C_block = malloc_matrix(N/(num_processes-1), N);
	if(rank==0)
	{
		A = malloc_matrix(N, M);
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

		printf("================== MATRIX A ===============\n");
		printMatrix(A, N, M);
		printf("================== MATRIX B ===============\n");
		printMatrix(B, M, N);
	}

	MPI_Bcast(B, N*M, MPI_FLOAT, 0, comm);
	MPI_Scatter(A, N*M/(num_processes), MPI_FLOAT, A_block, N*M/(num_processes), MPI_FLOAT, 0, comm);
	
	Matrix_Multiply(A_block, B, C_block, N/(num_processes), M, N);

	// printf("================MATRIX C_BLOCK from P%d============\n", rank);
	// printMatrix(C_block, N/(num_processes), N);


	if(rank==0) MPI_Gather(C_block, N*M/(num_processes), MPI_FLOAT, C, N*M/(num_processes), MPI_FLOAT, 0, comm);
	else {MPI_Gather(C_block, N*M/(num_processes), MPI_FLOAT, NULL, N*M/(num_processes), MPI_FLOAT, 0, comm);}
	
	if(rank==0) {
		float *D;
		D = malloc_matrix(N, N);
		for(int d = 0;d<(N*N);d++){
			D[d] = (float)0;
		}
		Matrix_Multiply(A, B, D, N, M, N);
		printf("================== MATRIX C_SERIAL ===============\n");
		printMatrix(D, M, N);

		printf("================== MATRIX C_PARALLEL ===============\n");
		printMatrix(C, M, N);
		
		printf("%d\n", isEqual(C, D, N*N));
	}

	MPI_Finalize();
	return 0;
}