all:
	/usr/local/bin/g++-9 -std=c++14 -O3 -fopenmp -o forth.out pthread.cpp

clean:
	rm *.out
