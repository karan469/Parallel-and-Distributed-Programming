#include "libs.h"

int main(int argc, char const *argv[])
{
	int num_processes, rank;
	MPI_Init(NULL, NULL);

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &num_processes);

	srand(time(NULL));
	double local_var = rand()/(double)RAND_MAX * 10;
	double num;
	MPI_Status status;
	if(rank == 0){
		MPI_Send(&local_var, 1, MPI_INT, 1, 1, comm);
		printf("Sended to 1: %.2f\n", local_var);
		
		MPI_Send(&local_var, 1, MPI_INT, 2, 2, comm);
		printf("Sended to 2: %.2f\n", local_var);

		MPI_Barrier(comm);

		MPI_Recv(&num, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
		printf("Master recieved from %d: %.2f\n",status.MPI_SOURCE, num);

		MPI_Recv(&num, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
		printf("Master recieved from %d: %.2f\n",status.MPI_SOURCE, num);
	} else if(rank == 1){
		MPI_Recv(&num, 1, MPI_INT, 0, MPI_ANY_TAG, comm, &status);
		printf("1 recieved from 0: %.2f\n", num);
		num = num * 2;
		MPI_Send(&num, 1, MPI_INT, 0, 90, comm);
	} else if(rank ==2){
		MPI_Recv(&num, 1, MPI_INT, 0, MPI_ANY_TAG, comm, &status);
		printf("2 recieved from 0: %.2f\n", num);
		num = num * 3;
		MPI_Send(&num, 1, MPI_INT, 0, 99, comm);
	}

	MPI_Finalize();

	return 0;
}