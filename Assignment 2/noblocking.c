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

	MPI_Status status;
    	MPI_Request request;


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

		// Printing Original matrices A and B
		// printf("================== MATRIX A ===============\n");
		// printMatrix(A, N, M);
		// printf("================== MATRIX B ===============\n");
		// printMatrix(B, M, N);
		// printMatrix(C, N, N);
		
		double start = MPI_Wtime();

		// Sending message to every other slave process their part of matrix
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

			MPI_Irsend(A_block, (int)(N*M/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*13, comm,&request);
			MPI_Irsend(B_block, (int)(M*N/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*97, comm,&request);

            		MPI_Wait(&request,&status);
		}

		C_block = malloc_matrix(N, N);
		for(int c = 0;c<(N*N);c++){
			C_block[c] = (float)0;
		}

		for(int f=0;f<num_processes-1;f++){
			MPI_Irecv(C_block, N*N, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &request);
            		MPI_Wait(&request,&status);
			addMatrices(C, C_block, C, N*N);
		}

		// Checking from Actual answer
		// float *D;
		// D = malloc_matrix(N, N);
		// for(int d = 0;d<(N*N);d++){
		// 	D[d] = (float)0;
		// }
		// Matrix_Multiply(A, B, D, N, M, N);
		// printf("%d\n", isEqual(C, D, N));   
        
        double end = MPI_Wtime();

        double duration = (float)end-start;

        printf("[NON-BLOCKING P2P] Time it took to run the process is %0.4fs\n",duration);
	}

	else if (rank>0)
	{
		// Allocating local block matrices for each process
		A_block = malloc_matrix(N, M/(num_processes-1));
		B_block = malloc_matrix(M/(num_processes-1), N);
		C_block = malloc_matrix(N,N);

		// Recieving with non-blocking mechanism
		MPI_Status status;
		MPI_Irecv(A_block, INT_MAX, MPI_FLOAT, 0, (rank)*13, comm, &request);
		MPI_Irecv(B_block, INT_MAX, MPI_FLOAT, 0, (rank)*97, comm, &request);

        	MPI_Wait(&request,&status);

		for(int c = 0;c<(N*N);c++){
			C_block[c] = (float)0;
		}

		// Matrix multiply serially on each process on sub matrices
		Matrix_Multiply(A_block, B_block, C_block, N, M/(num_processes-1), N);
		MPI_Isend(C_block, N*N, MPI_FLOAT, 0, (rank)*113, comm,&request);

			// Explicit wait before master process (Process 0) is ready
       	 	MPI_Wait(&request,&status);

	}

	MPI_Finalize(); // Ending the MPI program

	return 0;
}
