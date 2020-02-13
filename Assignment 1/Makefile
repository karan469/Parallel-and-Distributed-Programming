pthread:
	/usr/local/bin/g++-9 -std=c++11 -O3 -fopenmp -o pthread.out pthread.cpp
serial:
	/usr/local/bin/g++-9 -std=c++11 -O1 -fopenmp -o forth.out forth.cpp
openmp:
	/usr/local/bin/g++-9 -std=c++11 -O3 -fopenmp -o forth.out forth.cpp
all:
	/usr/local/bin/g++-9 -std=c++11 -O3 -fopenmp -o pthread.out pthread.cpp; g++ -std=c++11 -O3 -fopenmp -o forth.out forth.cpp
clean:
	rm *.out
