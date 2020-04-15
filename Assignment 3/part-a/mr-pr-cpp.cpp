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
#define ACC_TOLERANCE 0.0001
#define MAX_ITERATIONS 10000
#define alpha 0.5
using namespace std;

map<int, vector<int>> graph;
double pr[1000000];
int max_node_num=3;

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
    print("========================");
    for(int i=0;i<=max_node_num;i++){
        print(i<<": "<<pr[i]);
    }
    
}

string check_pr_status(double *old_pr){
    for(auto itr = graph.begin();itr!=graph.end();itr++){
        int node = itr->first;
        vector<int> outlinks = itr->second;
        double sum = 0;
        for(int i=0;i<outlinks.size();i++){
            if(outlinks[i]==-1) continue;
            auto it = graph.find(outlinks[i]);
            print(outlinks[i])
            assert(it!=graph.end());
            sum += old_pr[outlinks[i]]/(it->second.size());
        }
        sum = alpha*sum + (1-alpha)/(max_node_num+1);
        cout<<"[Node: "<<node<<"]"<<"Calculated now: "<<sum<<", In Table: "<<old_pr[node]<<endl;
        if(abs(sum-pr[node])>0.001){return "false";}
    }
    return "true";
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
        return (key * step_ <= last_);
    }

    bool const get_data(typename MapTask::key_type const &key, typename MapTask::value_type &value)
    {
        typename MapTask::value_type val;

        val.first  = first_ + (key * step_);
        val.second = std::min(val.first + step_ - 1, last_);
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
        // print("Input Value to map_task: ("<<value.first<<", "<<value.second<<")");
        for(int i=value.first;i<=value.second;i++){
            MAP_TYPE::iterator it = graph.find(i);
            if(it==graph.end()){
                print(i<<" ERROR");
                exit(1);
            }
            vector<int> l = it->second;
            
            if(l.size()==0) continue;
            
            for(int j=0;j<l.size();j++){
                
                if(l[j]==-1){break;}

                // if(l[j]>max_node_num) {print("HUSH "<<l[j]<<" "<<max_node_num); continue;}
                
                assert(l[j]<=max_node_num);
                
                // print("Output from map_task: Key="<<l[j]<<", Value="<<pr[j]/l.size());
                
                runtime.emit_intermediate(l[j], pr[j]/l.size());
            }        
        }
    }
};

struct reduce_task : public mapreduce::reduce_task<int, double>
{
    template<typename Runtime, typename It>
    void operator()(Runtime &runtime, key_type const &key, It it, It ite) const
    {
        value_type sum = 0;
        // cout<<"Input to reduce_task: Key="<<key<<" Value = [";
        for(auto itr=it;itr!=ite;itr++){
            // cout<<*itr;
            sum += *itr;
        }
        // print(sum);
        // cout<<"]\n";
        assert(key<=max_node_num);
        
        // assert(key!=0);
        pr[key] = (1-alpha) + alpha*(sum);
        // print("set for "<<key);
        // cout<<"Output from reduce_task: Key = "<<key<<", Value="<<double((1-alpha) + alpha*(sum))<<endl;
        
        runtime.emit(key, ((1-alpha)+ alpha*(sum)));
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
    max_node_num = 0;
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
    // File reading end

    // print("max_node_num: "<<max_node_num);

    for(int i=0;i<=max_node_num;i++){
        if(graph.find(i)==graph.end()) graph[i].push_back(-1);
    }

    // print_graph(graph);

    for(int i=0;i<=max_node_num;i++){
        pr[i] = (1);
    }

    // if (argc > 2)
    //     spec.map_tasks = std::max(1, atoi(argv[2]));
    spec.map_tasks = 1;
    int reduce_tasks = 1;
    // if (argc > 3)
    //     reduce_tasks = atoi(argv[3]);
    // else
    //     reduce_tasks = std::max(1U, std::thread::hardware_concurrency());
    spec.reduce_tasks = reduce_tasks;
    pagerank::job::datasource_type number_source(0, max_node_num, (max_node_num+1)/reduce_tasks);
    
    // std::cout <<"\nCalculating Page rank " << filename << " ..." <<std::endl;
    
    pagerank::job job(number_source, spec);
    mapreduce::results result;

    // main ops started
    int iter = 0;
    double sum_all_pr = 0;
    
    double old_pr[max_node_num+1];
    for(int i=0;i<=max_node_num;i++){
        old_pr[i] = 1;
    }
    double suma = 1e5;
    // print("HOLA");
    // print_pr();
    while((abs(suma)>ACC_TOLERANCE) && iter<MAX_ITERATIONS)
    // for(int i=0;i<12;i++)
    {
        iter += 1;
        #ifdef _DEBUG
            job.run<mapreduce::schedule_policy::sequential<pagerank::job> >(result);
        #else
            job.run<mapreduce::schedule_policy::cpu_parallel<pagerank::job> >(result);
        #endif

        for (auto it=job.begin_results(); it!=job.end_results(); ++it){
            // print(it->first);
            pr[it->first] = it->second;
        }

        // print_pr();
        
        suma = 0;
        for(int k=0;k<=max_node_num;k++){
            suma += (old_pr[k]-pr[k]);
        }

        for(int k=0;k<=max_node_num;k++){
            old_pr[k] = pr[k];
        }

        // cout<<check_pr_status(old_pr)<<endl;
    }

    // std::cout <<"\nMapReduce finished in " << result.job_runtime.count() << " with " << std::distance(job.begin_results(), job.end_results()) << " results" << std::endl;
    // print("Num of Iterations: "<<iter);
    
    double last_sum = 0;
    for(int i=0;i<=max_node_num;i++){
        last_sum += old_pr[i];
    }

    for(int i=0;i<=max_node_num;i++){
        print(i<<" = "<<pr[i]);
    }
    
    print("s = "<<last_sum);

	return 0;
}