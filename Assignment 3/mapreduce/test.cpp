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

map<int, vector<int>> graph;

void update_graph(vector<char> v){
	for(int i=0;i<v.size();i+=4){
		int a, b;
		a = int(v[i]-'0');
		b = int(v[i+2]-'0');
		graph[a].push_back(b);
	}
}

void print_graph(map<int, vector<int> > graph){
	for(auto itr = graph.begin(); itr!=graph.end(); itr++){
		cout<<itr->first<<": [";
		for(int i=0;i<itr->second.size();i++){
			cout<<itr->second[i]<<", ";
		}
		cout<<"]"<<endl;
	}
}

int main(int argc, char const *argv[])
{
	// Reading file
	std::vector<char> v;
	if (FILE *fp = fopen("diamond.txt", "r"))
	{
		char buf[1024];
		while (size_t len = fread(buf, 1, sizeof(buf), fp))
			v.insert(v.end(), buf, buf + len);
		fclose(fp);
	}
	update_graph(v);
	print_graph(graph);

	// graph filled



	return 0;
}