all:
	g++ -std=c++14 -O3 -fopenmp -o forth.out forth.cpp
clean:
	rm *.out