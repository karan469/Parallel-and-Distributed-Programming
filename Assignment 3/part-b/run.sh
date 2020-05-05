mpic++ -o mr-pr-mpi.o  mr-pr-mpi.cpp; mpirun -n 3 ./mr-pr-mpi.o $1.txt -o $1-pr-mpi.txt
