#include "libs.h"

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

	// if (num_processes < 2) {
	//     fprintf(stderr, "Must use atleast two processes for this example\n");
	//     MPI_Abort(comm, 1);
	// }

	if(N%num_processes!=0){
		fprintf(stderr, "Cant divide work properly to all processes\n");
	    MPI_Abort(comm, 1);
	}


	B = malloc_matrix(M, N);
	A_block = malloc_matrix(N/(num_processes), M);
	C_block = malloc_matrix(N/(num_processes), N);

	//Enter master process
	if(rank==0)
	{
		A = malloc_matrix(N, M);
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
	}

	double start = MPI_Wtime();

	//Broadcasts/Sends Matrix B to all processes
	MPI_Bcast(B, N*M, MPI_FLOAT, 0, comm);

	//Scatters/Divides A into blocks sending them individually to all processes
	MPI_Scatter(A, N*M/(num_processes), MPI_FLOAT, A_block, N*M/(num_processes), MPI_FLOAT, 0, comm);
	
	Matrix_Multiply(A_block, B, C_block, N/(num_processes), M, N);

	// printf("================MATRIX C_BLOCK from P%d============\n", rank);
	// printMatrix(C_block, N/(num_processes), N);

	//Gathers computations and concatenate answers with each other
	if(rank==0) MPI_Gather(C_block, N*N/(num_processes), MPI_FLOAT, C, N*N/(num_processes), MPI_FLOAT, 0, comm); // if master, concatenate computed block with output matrix
	else {MPI_Gather(C_block, N*N/(num_processes), MPI_FLOAT, NULL, N*N/(num_processes), MPI_FLOAT, 0, comm);}  // if not master, receive only the block
	
	// Do wee need this?
	// MPI_Barrier(comm);

	if(rank==0){
		double end = MPI_Wtime();

	    double duration = (float)end-start;

	    printf("[COLLECTIVE] Time it took to run the process is %0.4fs\n",duration);
	}

	// Checking answer with serial
	if(rank==0) {
		float *D;
		D = malloc_matrix(N, N);
		for(int d = 0;d<(N*N);d++){
			D[d] = (float)0;
		}
		Matrix_Multiply(A, B, D, N, M, N);
		// printf("================== MATRIX C_SERIAL ===============\n");
		// printMatrix(D, N, N);

		// printf("================== MATRIX C_PARALLEL ===============\n");
		// printMatrix(C, N, N);
		
		printf("%d\n", isEqual(C, D, N*N));
	}

	//Ends MPI Environment
	MPI_Finalize();
	return 0;
}