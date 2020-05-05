// MapReduce word frequency example in C++
// Syntax: wordfreq file1 dir1 file2 dir2 ...
// (1) read all files and files in dirs
// (2) parse into words separated by whitespace
// (3) count occurrence of each word in all files
// (4) print top 10 words

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include "stdc++.h"

void input(int, MAPREDUCE_NS::KeyValue *, void *);
void mymap(int, MAPREDUCE_NS::KeyValue *, void *);
void myreduce(char *, int, char *, int, int *, MAPREDUCE_NS::KeyValue *, void *);
void change_I(uint64_t, char *, int, char *, int, MAPREDUCE_NS::KeyValue *, void *);

struct PAGE_RANK{
	std::string filename;
	std::vector<std::vector<std::pair<unsigned, double> > > H;
	unsigned n;
	double* I;
	std::vector<unsigned> A;
	int nprocs;
	double alpha;
	double diff;
	double avalue;
	double one;
};

int main(int narg, char **args)
{
	MPI_Init(&narg,&args);

	int me,nprocs, loop_count=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&me);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	if(narg <= 2){
		if(me == 0){
			printf("Syntax: wordfreq file1 file2\n");
		}
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	PAGE_RANK page_rank;
	page_rank.filename = args[1];
	page_rank.nprocs = nprocs;
	std::string file2 = args[2];
	// MPI_Barrier(MPI_COMM_WORLD);
	// double tstart = MPI_Wtime();


	
	MAPREDUCE_NS::MapReduce *mr = new MAPREDUCE_NS::MapReduce(MPI_COMM_WORLD);
	// mr->timer = 1;

	mr->map(nprocs, &input, &page_rank); // 1 or nprocs?
	delete mr;


	MPI_Barrier(MPI_COMM_WORLD);
	double tstart = MPI_Wtime();
L:
	loop_count ++;
	MPI_Barrier(MPI_COMM_WORLD);
	mr = new MAPREDUCE_NS::MapReduce(MPI_COMM_WORLD);
	// mr->timer = 1;
	page_rank.diff = 0;
	page_rank.avalue = 0;
	page_rank.one = 0;
	page_rank.alpha = 0.85;
	mr->map(nprocs, &mymap, &page_rank);
	mr->convert();
	// mr->collate(NULL);
	mr->reduce(myreduce, &page_rank);
	mr->gather(1);
	// int key=0;
	// mr->collapse((char *)&key, sizeof(int));
	mr->map(mr, change_I, &page_rank);
	delete mr;
	MPI_Bcast(&page_rank.diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(page_rank.I, page_rank.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(page_rank.diff>0.0000001)
		goto L;

	MPI_Barrier(MPI_COMM_WORLD);
	double tstop = MPI_Wtime();
	if(me==0){
		std::fstream fout(file2, std::ios::out);
		for(unsigned i=0; i<page_rank.n; ++i){	
			fout<<i<<" = "<<page_rank.I[i]<<"\n";
		}
		std::cout<<"\n"<<"TIME TAKEN in " <<loop_count<<" loops : "<<tstop-tstart<<"\n\n";
	}
	MPI_Finalize();
}

void mymap(int itask, MAPREDUCE_NS::KeyValue *kv, void *ptr){
	PAGE_RANK *P = (PAGE_RANK *) ptr;
	int n = P->n;
	int nprocs = P->nprocs;
	int start=itask*(n/nprocs), end;
	if(itask==P->nprocs-1){
		end = n;
	}
	else{
		end = (itask+1)*(n/nprocs);
	}
	for(unsigned i=start; i<end; ++i){
		kv->add((char *) &i, sizeof(unsigned), (char *)&(P->H[i]), P->H[i].size()*sizeof(std::pair<unsigned, double>));
	}
	if(itask == 0){
		for(unsigned i=0; i<P->n; ++i){	
	    	P->one += P->I[i]/((double)P->n);
	    }
	    for(unsigned i=0; i<P->A.size(); ++i){
	    	P->avalue += P->I[P->A[i]]/((double)P->n);
	    }
	}
}

void myreduce(char *key, int keybytes, char *multivalue,
	 int nvalues, int *valuebytes, MAPREDUCE_NS::KeyValue *kv, void *ptr)
{
	double res = 0;
	if(*valuebytes == 0){
		kv->add(key, keybytes, (char *)&res, sizeof(double));
		return;
	}
	PAGE_RANK *P = (PAGE_RANK *) ptr;
	std::vector<std::pair<unsigned, double> > *rowptr = (std::vector<std::pair<unsigned, double> > *) multivalue;
	// std::cout<<*((unsigned *)key)<<": "<<*valuebytes<<"\n";
	std::vector<std::pair<unsigned, double> > val = *rowptr;
	for(unsigned i=0; i<val.size(); ++i){
    	unsigned c = val[i].first;
    	double v = val[i].second;
    	res += P->I[c]*v;
	}
	kv->add(key, keybytes, (char *)&res, sizeof(double));
}

void change_I(uint64_t itask, char *key, int keybytes, char *value,
		int valuebytes, MAPREDUCE_NS::KeyValue *kv, void *ptr)
{
	PAGE_RANK *P = (PAGE_RANK *) ptr;
	double *res = (double *) value;
	unsigned *index = (unsigned *) key;
	double a = P->alpha;
	P->diff += std::abs(a*(*res) + a*P->avalue + (1-a)*P->one - P->I[*index]);
	P->I[*index] = a*(*res) + a*P->avalue + (1-a)*P->one;
	// unsigned *index = (unsigned *) value;
	// double *res = (double *)(index+1);
	// unsigned i=0;
	// while(i<valuebytes){
	// 	unsigned *index = (unsigned *) (value+i);
	// 	double *res = (double *)(index+1);
	// 	i += sizeof(index)+sizeof(res);
	// 	P->diff += std::abs(a*(*res) + a*P->avalue + (1-a)*P->one - P->I[*index]);
	// 	P->I[*index] = a*(*res) + a*P->avalue + (1-a)*P->one;
	// }
}

void input(int itask, MAPREDUCE_NS::KeyValue *kv, void *ptr){
	PAGE_RANK *P = (PAGE_RANK *) ptr;
	std::fstream fin;
	fin.open(P->filename, std::ios::in);
	P->n=0;
	while(!fin.eof()){
		unsigned x;
		fin>>x;
		P->n = std::max(P->n, x);
	}
	P->n ++;
	std::vector<double> colSums(P->n);
	colSums.assign(P->n, 0);
	fin.close();
	fin.open(P->filename, std::ios::in);
	std::vector<std::pair<unsigned, double> > H_dum;
	P->H.assign(P->n, H_dum);
	unsigned i, j;
	while(fin>>i){
		fin>>j;
		P->H[j].push_back(std::make_pair(i, 1.0));
		colSums[i] += 1;
	}
	for(unsigned i=0; i<P->n; ++i){
		if(colSums[i]==0){
			P->A.push_back(i);
		}
		for(unsigned j=0; j<P->H[i].size(); ++j){
			std::pair<unsigned, double> val = P->H[i][j];
			unsigned col = val.first;
			double change = val.second/colSums[col];
			P->H[i][j] = std::make_pair(col, change);
		}
	}
	fin.close();
	P->I=new double[P->n];
	for(unsigned i=0; i<P->n; ++i)
		P->I[i] = 0;
	P->I[0] = 1;
}

