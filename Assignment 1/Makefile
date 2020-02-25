pthread:
	g++ -std=c++11 -O3 -fopenmp -o pthread.out pthread.cpp
serial:
	g++ -std=c++11 -O1 -fopenmp -o forth.out forth.cpp
openmp:
	g++ -std=c++11 -O3 -fopenmp -o forth.out forth.cpp
all:
	g++ -std=c++11 -O3 -fopenmp -o pthread.out pthread.cpp; g++ -std=c++11 -O3 -fopenmp -o forth.out forth.cpp
clean:
	rm *.out
