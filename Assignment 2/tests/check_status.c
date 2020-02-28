#include "libs.h"

int main(int argc, char const *argv[])
{
	const int MAX_NUM = 100;
	int numbers[MAX_NUM];
	for(int i=0;i<MAX_NUM;i++){
		numbers[i] = i+1;
	}
	/* code */
	int rank, num_processes;
	MPI_Init(NULL, NULL);

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &num_processes);

	if (num_processes < 2) {
	    fprintf(stderr, "Must use atleast two processes for this example\n");
	    MPI_Abort(comm, 1);
	}


	// printf("Process %d out of %d\n", rank, num_processes);

	if(rank == 0){
		srand(time(NULL));
		int amount = (rand()/(float)RAND_MAX) * MAX_NUM;

		MPI_Send(numbers, amount, MPI_INT, 1, 99, comm);
		printf("> 0 sent %d numbers to Process 1\n", amount);
	} else if(rank == 1){
		int numbers_recv[MAX_NUM];
		MPI_Status status;
		MPI_Recv(numbers_recv, MAX_NUM, MPI_INT, 0, 99, comm, &status);

		int num_amount;

		MPI_Get_count(&status, MPI_INT, &num_amount);
		
		for(int i = 0;i<num_amount; i++){
			printf("%d ", numbers_recv[i]);
		}
		printf("\n> 1 recieved %d numbers from process %d with message id %d\n", num_amount, status.MPI_SOURCE, status.MPI_TAG);
	}

	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;
}