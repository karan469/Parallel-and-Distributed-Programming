#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv){
	int np, mype, ierr;

	MPI_Init(&argc, &argv);
	// printf("%d\n", ierr);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	// printf("%d\n", ierr);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	// printf("%d\n", ierr);

	printf("total processes = %d | my rank = %d\n", np, mype);

	MPI_Finalize();
}