Description
==

This assignment contains 3 parts in which we are calculating pagerank of webpages using mapreduce libraries written in c++ using pthreads and mpi.

Each folder contains script.py to compile and run code and calculate time of execution for each benchmark given in ./test/ folder in each part.

```
	python3 script.py > exec_time_report.csv
```

For each input file, ./run.sh is run with 4 processes and one arguement.
```
	./run.sh {filename-without-extension}
```

run.sh script compiles and runs the c++ program with 2 arguements: input filename and output filename.
```
	./mr-pr-cpp.o ${filename}.txt -o ${filename}-pr-cpp.txt // For Part-A
	./mr-pr-mpi.o ${filename}.txt -o ${filename}-pr-mpi.txt // For Part-B
```
