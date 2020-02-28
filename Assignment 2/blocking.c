#include "libs.h"
#define BUFSIZE INT_MAX

//Input/Output Matrices
float *A, *B, *C;
float *A_block, *B_block, *C_block;

int main(int argc, char const *argv[])
{
	const int M = 32;
	
	int N =  atoi(argv[1]);
	

	int rank, num_processes;

	//Initialises Open MPI environment
	MPI_Init(NULL, NULL);

	//Designate number of process and rank to each process
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &num_processes);

	if (num_processes < 2) {
	    fprintf(stderr, "Must use atleast two processes for this example\n");
	    MPI_Abort(comm, 1);
	}

	//Enter master process
	if(rank == 0)
	{
		int a_mssg_id = 0;
		int b_mssg_id = 0;

		A = malloc_matrix(N, M);
		B = malloc_matrix(M, N);
		C = malloc_matrix(N, N);

		//Random Initialisation matrices to multiply
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
		

		//Dividing Input Matrices into blocks to send to slave processes
		for(int f=0;f<num_processes-1;f++){
			A_block = malloc_matrix(N, M/(num_processes-1));
			B_block = malloc_matrix(M/(num_processes-1), N);

			for(int i=0;i<N*M/(num_processes-1);i++){
				B_block[i] = B[i + (N*M*f)/(num_processes-1)];
			}

			// Matrix_Multiply(A, B, C, N, M, N);
			// printf("REAL ANS process  %d\n", f+1);
			// printMatrix(C, N, N);
			// printf("ENDS %d\n", f+1);

			int counter = 0;
			for(int i=0;i<N;i++){
				for(int j = i*M + f*(M/(num_processes-1)); j<i*M + (int)(M/(num_processes-1)) + f*(M/(num_processes-1));j++){
					A_block[counter++] = A[j];
				}
			}

			//Blocking sending of blocks of matrices to slave
			MPI_Rsend(A_block, (int)(N*M/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*13, comm);
			MPI_Rsend(B_block, (int)(M*N/(num_processes-1)), MPI_FLOAT, f+1, (f+1)*97, comm);
		}

		C_block = malloc_matrix(N, N);
		for(int c = 0;c<(N*N);c++){
			C_block[c] = (float)0;
		}

		MPI_Status status;
		for(int f=0;f<num_processes-1;f++){
			//Receives multiplied blocks from slave processes
			MPI_Recv(C_block, N*N, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
			addMatrices(C, C_block, C, N*N);
		}
		// printf("=======MPI Ans=========\n");
		// printMatrix(C, N, N);
		float *D;
		D = malloc_matrix(N, N);
		for(int d = 0;d<(N*N);d++){
			D[d] = (float)0;
		}
		// printf("=======ACTUAL Ans=========\n");
		Matrix_Multiply(A, B, D, N, M, N);
		// printMatrix(D, N, N);

		//Compare serial and parallel answer
		printf("%d\n", isEqual(C, D, N));
	}
	else if (rank>0)
	{
		// printf("process %d started..\n", rank);
		A_block = malloc_matrix(N, M/(num_processes-1));
		B_block = malloc_matrix(M/(num_processes-1), N);
		C_block = malloc_matrix(N,N);

		MPI_Status status;

		//Receive small blocks of input matrices
		MPI_Recv(A_block, INT_MAX, MPI_FLOAT, 0, (rank)*13, comm, &status);
		MPI_Recv(B_block, INT_MAX, MPI_FLOAT, 0, (rank)*97, comm, &status);

		for(int c = 0;c<(N*N);c++){
			C_block[c] = (float)0;
		}

		//Multiply small blocks 
		Matrix_Multiply(A_block, B_block, C_block, N, M/(num_processes-1), N);
		// printMatrix(C_block, N, N);

		//Send the computation to master process
		MPI_Send(C_block, N*N, MPI_FLOAT, 0, (rank)*113, comm);
	}

	//Ends the MPI environment
	MPI_Finalize();
	return 0;
}
