compile:
	mpicc -o $(input).o $(input).c
run:
	mpirun -n $(num_proc) --use-hwthread-cpus $(input).o $(N)

