#include "mpi.h"
#include "stdc++.h"
#include<unordered_map>
#define BUFSIZE INT_MAX
#define comm MPI_COMM_WORLD
#define print(x) cout<<x<<"\n";
using namespace std::chrono;
using namespace std;


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

void print_graph(unordered_map<int,vector<int> >& graph){
    for(auto it=graph.begin();it!=graph.end();it++){
        cout << it->first << " : ";
        vector<int> graph_vec = it->second;
        int sizes = it->second.size();
        for(int i=0;i<sizes;i++){
            cout << graph_vec[i] << " ";
        }
        cout << endl;
    }
}

multimap<float, int> invert(unordered_map<int, float>& mymap)
{
	multimap<float, int> multiMap;

	unordered_map<int, float> :: iterator it;
	for (it=mymap.begin(); it!=mymap.end(); it++)
	{
		multiMap.insert(make_pair(it->second, it->first));
	}

	return multiMap;
}

unordered_map<int,float> pr_init(unordered_map<int,vector<int> >& graph){
    unordered_map<int,float> pr_map;

    for(auto it=graph.begin();it!=graph.end();it++){
        pr_map[it->first] = (1/(float)graph.size());
    }
    return pr_map;
}

void mpi_pr(unordered_map<int,vector<int> >& graph, unordered_map<int,float>& pr_map,int rank,int num_processes){

    unordered_map<int,vector<float> > map_input;
    int id;

    if(rank==0){
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

        for(auto it=graph.begin();it!=graph.end();it++){
            id = it->first;
            vector<float> in_neg = map_input[id];
            MPI_Rsend(&id, 1, MPI_INT, ((it->first)%(num_processes-1))+1, ((it->first+1)*10),comm);
            MPI_Rsend(&in_neg[0], in_neg.size(), MPI_FLOAT, ((it->first)%(num_processes-1))+1, ((it->first+1)*10)+1, comm);
        }

        MPI_Status status;

        for(auto it=graph.begin();it!=graph.end();it++){
            float pr=0;
		    MPI_Recv(&pr, 1, MPI_FLOAT, ((it->first)%(num_processes-1))+1, (it->first)+1, comm, &status);
            // cout << "Message " << pr << "  Received by master from id " << it->first << endl;
            pr_map[it->first] = pr;
        }
    }

    else if(rank>0){
        MPI_Status status;
        for(int i=rank;i<=graph.size();i+=(num_processes-1)){
            pair<int,vector<float> > pi;
            vector<float> neg_values;
		    MPI_Recv(&id, INT_MAX, MPI_INT, 0, i*10, comm, &status);
            neg_values.resize(graph.size());
		    MPI_Recv(&neg_values[0], INT_MAX, MPI_FLOAT, 0, i*10+1, comm, &status);
            float sum = 0;
            for(auto it=neg_values.begin();it!=neg_values.end();it++){
                sum += *it;
            }

            float temp_pr = (float)(0.15/graph.size()) + (float)(0.85*sum);
            MPI_Send(&temp_pr, 1, MPI_FLOAT, 0, i, comm);
        }
    }
}

int main(int argc, char *argv[]){
  int num_processes = 5;
  int rank;
  string filename = "./input.txt";
  if(argc>1){
    filename = argv[1];
  }
  unordered_map<int,vector<int> > graph = read_graph(filename);
  MPI_Init(NULL, NULL);

	//Designate number of process and rank to each process
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &num_processes);

	if (num_processes < 1) {
	    cerr << "Must use atleast two processes for this example\n";
	    MPI_Abort(comm, 1);
	}

    unordered_map<int,float> pr_map = pr_init(graph);

    auto start = high_resolution_clock::now();
    for(int i=0;i<60;i++){
        mpi_pr(graph,pr_map,rank,num_processes);
    }

    if(rank==0){
        multimap<float, int> inverted = invert(pr_map);
        float sum = 0;
        vector<pair<int, float> > finalList;
        for(auto it=inverted.begin();it!=inverted.end();it++){
            finalList.push_back(make_pair(it->second, it->first));
            sum += it->first;
        }
        sort(finalList.begin(), finalList.end());

        ofstream fout;
        string outfilename = "file-pr-mpi.txt";
        if(argc>3){
          outfilename = argv[3];
        }
        fout.open(outfilename);

        for(int i=0;i<finalList.size();i++){
          fout<<finalList[i].first<<" = "<<finalList[i].second/(float)sum<<'\n';
          // cout<<finalList[i].first<<" = "<<finalList[i].second/(float)sum<<'\n';
        }
        // cout<<"s = "<<1<<'\n';
        fout<<"s = "<<1<<'\n';
        fout.close();
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        // print("** Time of execution = "<<duration.count()/(double)1e6<<" **");
        print(filename<<", "<<duration.count()/(double)1e6);
    }

	//Ends the MPI environment
	MPI_Finalize();


  return 0;
}
