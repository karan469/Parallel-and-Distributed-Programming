// Copyright (c) 2009-2016 Craig Henderson
// https://github.com/cdmh/mapreduce

// The prime number code is based on work from Christian Henning [chhenning@gmail.com]

#include <boost/config.hpp>
#if defined(BOOST_MSVC)
#   pragma warning(disable: 4127)

// turn off checked iterators to avoid performance hit
#   if !defined(__SGI_STL_PORT)  &&  !defined(_DEBUG)
#       define _SECURE_SCL 0
#       define _HAS_ITERATOR_DEBUGGING 0
#   endif
#endif

#include "./include/mapreduce.hpp"
#include <iostream>
#define print(x) std::cout<<x<<endl;
#define MAP_TYPE map<int, vector<int>>
#define ACC_TOLERANCE 0.005
using namespace std;

map<int, vector<int>> graph;
float pr[1000];
int max_node_num = 11;

void print_graph(map<int, vector<int> > graph){
    for(auto itr = graph.begin(); itr!=graph.end(); itr++){
        cout<<itr->first<<": [";
        for(int i=0;i<itr->second.size();i++){
            cout<<itr->second[i]<<", ";
        }
        cout<<"]"<<endl;
    }
}

void print_pr(){
    for(int i=0;i<max_node_num;i++){
        print(pr[i]);
    }
}

namespace pagerank {

template<typename MapTask>
class number_source : mapreduce::detail::noncopyable
{
  public:
    number_source(long first, long last, long step)
      : sequence_(0), first_(first), last_(last), step_(step)
    {

    }

    bool const setup_key(typename MapTask::key_type &key)
    {
        key = sequence_++;
        // print("setup_key called "<<key<<" x "<<step_<<" x "<<last_<<" "<<"should proceed? "<<(key * step_ <= last_));
        return (key * step_ <= last_);
    }

    bool const get_data(typename MapTask::key_type const &key, typename MapTask::value_type &value)
    {
        typename MapTask::value_type val;

        val.first  = first_ + (key * step_);
        val.second = std::min(val.first + step_ - 1, last_);
        // cerr<<"VAL1: ["<<val.first<<", "<<val.second<<"]"<<endl;
        std::swap(val, value);
        return true;
    }

  private:
    long       sequence_;
    long const step_;
    long const last_;
    long const first_;
};

struct map_task : public mapreduce::map_task<int, pair<int, int> >
{
    template<typename Runtime>
    void operator()(Runtime &runtime, key_type const /*&key*/, value_type const &value) const
    {
        // print("VAL2: ["<<value.first<<" "<<value.second<<"]");
        for(int i=value.first;i<=value.second;i++){
            MAP_TYPE::iterator it = graph.find(i);
            if(it==graph.end()){
                print(i<<" ERROR");
                exit(1);
            }
            vector<int> l = it->second;
            for(int j=0;j<l.size();j++){
                if(l[j]==-1){continue;}
                // print("EMIT: {   "<<l[j]<<" | "<<pr[i]/l.size()<<"   }");
                runtime.emit_intermediate(l[j], pr[j]/l.size());
            }        
        }
    }
};

struct reduce_task : public mapreduce::reduce_task<int, float>
{
    template<typename Runtime, typename It>
    void operator()(Runtime &runtime, key_type const &key, It it, It ite) const
    {
        value_type sum = 0;
        // print("In reduce task: "<<key);
        for(auto itr=it;itr!=ite;itr++){
            // print(*itr);
            sum += *itr;
        }
        pr[key] = 0.15/max_node_num + 0.85*(sum);
        runtime.emit(key, (0.15 + 0.85*(sum)));
    }
};

typedef
mapreduce::job<pagerank::map_task,
               pagerank::reduce_task,
               mapreduce::null_combiner,
               pagerank::number_source<pagerank::map_task>
> job;

} // namespace pagerank

int main(int argc, char *argv[])
{
    

    mapreduce::specification spec;

    string filename = "diamond.txt";
    if (argc > 1)
        filename = string(argv[1]);

    // Reading file
    ifstream fin;
    fin.open(filename);
    int max_node_num = 0;
    int maxs = 0;
    if(fin){
        int s;
        while(fin>>s){
            int b;
            fin>>b;
            graph[s].push_back(b);
            if(s>maxs){
                maxs = s;
            }
            if(max_node_num<b){
                max_node_num = b;
            }
            if(max_node_num<s){
                max_node_num = s;
            }
        }
    }
    for(int i=0;i<=max_node_num;i++){
        if(graph.find(i)==graph.end()) graph[i].push_back(-1);
    }

    // print_graph(graph);
    // print(max_node_num);

    // print("PAGERANK BEFORE OPERATIONs-->");

    for(int i=0;i<=max_node_num;i++){
        pr[i] = 0.15/max_node_num;
    }
    
    // for(int i=0;i<=max_node_num;i++){
    //     print(i<<": "<<pr[i]);
    // }

    if (argc > 2)
        spec.map_tasks = std::max(1, atoi(argv[2]));

    int reduce_tasks;
    if (argc > 3)
        reduce_tasks = atoi(argv[3]);
    else
        reduce_tasks = std::max(1U, std::thread::hardware_concurrency());
    spec.reduce_tasks = reduce_tasks;

    // print("SPREC: "<<max_node_num<<" "<<max_node_num/reduce_tasks);
    pagerank::job::datasource_type number_source(0, max_node_num, (max_node_num+1)/reduce_tasks);

    // std::cout <<"\nCalculating Page rank " << filename << " ..." <<std::endl;
    
    // 
    pagerank::job job(number_source, spec);
    mapreduce::results result;


    float last_sum = 0;
    int iter = 0;
    for(int i=0;i<atoi(argv[4]);i++)
    // while(true)
    {
        iter += 1;
#ifdef _DEBUG
    job.run<mapreduce::schedule_policy::sequential<pagerank::job> >(result);
#else
    job.run<mapreduce::schedule_policy::cpu_parallel<pagerank::job> >(result);
#endif

    float sum_all_pr = 0;
    for(int k=0;k<=max_node_num;k++){
        sum_all_pr += pr[k];
    }

    // if((sum_all_pr-last_sum)*(sum_all_pr-last_sum) < ACC_TOLERANCE){break;}
    last_sum = sum_all_pr;
}
    // std::cout <<"\nMapReduce finished in " << result.job_runtime.count() << " with " << std::distance(job.begin_results(), job.end_results()) << " results" << std::endl;
    // print("Num of Iterations: "<<iter);
    for(int i=0;i<=max_node_num;i++){
        print(i<<" = "<<pr[i]/last_sum);
    }
    float temp = 0;
    for(int k=0;k<=max_node_num;k++){
        temp += pr[k]/last_sum;
    }
    print("s = "<<temp);


	return 0;
}

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
