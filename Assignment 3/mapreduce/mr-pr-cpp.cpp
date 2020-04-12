#include <bits/stdc++.h>
#include "./include/mapreduce.hpp"
#include <boost/config.hpp>
#include <fstream>
#include<stdlib.h>
#include <stdio.h>
#define print(x) std::cout<<x<<'\n'
#define min(x, y) -1?x<y:1
#define DAMPING_FACTOR 0.85

using namespace std;

map<int, float> pr;
vector<pair<int, vector<int> > > graph;

// void calc_pagerank()
void init_pr(int start, int end){
	pr.clear();
	for(int i=start;i<=end;i++){
		pr.insert({i, 0.18});
	}
}

void init_graph(int start, int end){
	srand(time(0));
	for(int i=start;i<=end;i++){
		int a = start + ( std::rand() % ( end/2 - start + 1 ) );
		int b = end/2 + ( std::rand() % ( end - end/2 + 1 ) );
		vector<int> outlinks;
		for(int j=a;j<=b;j++){
			outlinks.push_back(j);
		}
		graph.push_back({i, outlinks});
	}
}

void print_int_vector(vector<int> vec){
	for(int i=0;i<vec.size();i++){
		cout<<vec[i]<<" ";
	}
	cout<<endl;
}

void print_pr(){
	print("--- Printing Page Ranks ---");
	for(auto itr = pr.begin(); itr!=pr.end(); itr++){
		cout<<itr->first<<": "<<itr->second<<endl;
	}
}

void print_graph(){
	print("--- Printing Graphs ---");
	for(int i=0;i<graph.size();i++){
		cout<<graph[i].first<<": "<<endl;
		print_int_vector(graph[i].second);
	}
}

void pr_calculate(){
	for(int i=0;i<graph.size();i++){
		vector<int> outlinks = graph.second;
		int node_num = graph.first;
		auto map<int, float>::iterator itr = pr.find(node_num);
		float sum = (1-DAMPING_FACTOR)
	}
}

int main(int argc, char const *argv[])
{
	init_pr(1, 10);
	print_pr();
	init_graph(1, 10);
	print_graph();
	print("Hello world");
	
	return 0;
}