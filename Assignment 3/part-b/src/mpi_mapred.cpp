#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include "stdc++.h"
#include<unordered_map>
#define BUFSIZE INT_MAX
#define comm MPI_COMM_WORLD

// using namespace MAPREDUCE_NS;
using namespace std;

struct PageRank{
    unordered_map<int,float> page_rank;
    unordered_map<int,vector<int> > graph;
    unordered_map<int,vector<float> > map_input;
    int num_of_processes;
};

unordered_map<int,float> pr_init(unordered_map<int,vector<int> >& graph){
    unordered_map<int,float> pr_map;

    for(auto it=graph.begin();it!=graph.end();it++){
        pr_map[it->first] = (1/(float)graph.size());
    }

    return pr_map;
}

void mymap(int itask, MAPREDUCE_NS::KeyValue *kv, void *ptr){

    PageRank *p = (PageRank *) ptr;

    unordered_map<int,float> page_rank = p->page_rank;
    unordered_map<int,vector<int> > graph = p->graph;
    unordered_map<int,vector<float> > map_input = p->map_input;
    int num_of_processes = p->num_of_processes;

    for(auto it=graph.begin();it!=graph.end();it++){
        if((it->first)%(num_of_processes)==(itask)){       //if Node belongs to processor's part of graph
            kv->add((char *) &(it->first), sizeof(int), (char *) &((map_input[it->first])[0]), sizeof(float)*map_input[it->first].size());
        }
    }

}

void myreduce(char *key, int keybytes, char *multivalue,
	int nvalues, int *valuebytes, MAPREDUCE_NS::KeyValue *kv, void *ptr){

    int *index = (int *) key;
    vector<float> *inlinks = (vector<float> *) multivalue;
    PageRank *p = (PageRank *) ptr;

    float sum = 0;

    for(int i=0;i<(*inlinks).size();i++){
        sum += (*inlinks)[i];
    }

    float temp_pr = (float)(0.15/(float)(p->graph).size()) + (float)(0.85*sum);

    p->page_rank[*index] = temp_pr;

}


unordered_map<int,vector<int> > read_graph(string file){
    unordered_map<int,vector<int> > graph;

    string line;
    ifstream myfile(file);

    if(myfile.is_open()){
        while(getline(myfile,line)){
            istringstream ss(line);

            string word1;
            string word2;

            ss >> word1;
            ss >> word2;

            int id1,id2;

            id1 = stoi(word1);  
            id2 = stoi(word2);

            if(graph.find(id1)==graph.end()){
                vector<int> neg;
                neg.push_back(id2);
                graph[id1] = neg;
                if(graph.find(id2)==graph.end()){
                    vector<int> rand;
                    graph[id2] = rand;
                }
            }           
            else{
                graph[id1].push_back(id2);
                if(graph.find(id2)==graph.end()){
                    vector<int> rand;
                    graph[id2] = rand;
                }
            }

        }
    }

    return graph;
}

unordered_map<int,vector<float> > reverse_graph(unordered_map<int,vector<int> > &graph, unordered_map<int,float> &pr_map){
    unordered_map<int,vector<float> > map_input;

    for(auto it1=graph.begin();it1!=graph.end();it1++){
        vector<float> inlinks;
        for(auto it2=graph.begin();it2!=graph.end();it2++){
            vector<int> neg = it2->second;
            if(std::find(neg.begin(),neg.end(),it1->first)!=neg.end()){
                inlinks.push_back(it2->first);
            }
        }
        for(int i=0;i<inlinks.size();i++){
            inlinks[i] = ((float)pr_map[inlinks[i]]/(float)graph[inlinks[i]].size());
        }
        if (inlinks.size()==0) inlinks.push_back(0);
        map_input[it1->first] = inlinks;
    }

    return map_input;
}

int main(int narg, char **args)
{      
    MPI_Init(&narg,&args);

    int me,nprocs;
    nprocs = 5;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if (narg < 1) {
        printf("Syntax: {input file}\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();

    PageRank page_rank;
    page_rank.num_of_processes = nprocs;
    page_rank.graph = read_graph("./input.txt");
    page_rank.page_rank = pr_init(page_rank.graph);
    page_rank.map_input = reverse_graph(page_rank.graph,page_rank.page_rank);

    for(int i=0;i<10;i++){
        MPI_Barrier(MPI_COMM_WORLD);
        MAPREDUCE_NS::MapReduce *mr = new MAPREDUCE_NS::MapReduce(MPI_COMM_WORLD);

    	mr->map(nprocs, &mymap, &page_rank);
        mr->collate(NULL);
        mr->reduce(&myreduce, &page_rank);
        mr->gather(1);
        delete mr;
    }

    double time_stop = MPI_Wtime();

	MPI_Finalize();

    for(auto it=(page_rank.page_rank).begin();it!=(page_rank.page_rank).end();it++){
        cout << it->first << " : " << it->second << endl;
    }    

    cout << "Time taken is " << (time_stop-time_start) << endl;

    return 0;
}



