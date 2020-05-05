#include <boost/config.hpp>
#if defined(BOOST_MSVC)
#   pragma warning(disable: 4127)

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
#define MAX_ITERATIONS 1000
#define alpha 0.85
using namespace std::chrono;
using namespace std;

map<int, vector<int>> graph;
vector<double> pr(10000, 0.0f);
// vector<vector<double> > mm( 10000 , vector<double> (10000,0));
double mm[10000][10000];
vector<double> dangling(1000000, 0.0f);
double de;
int max_node_num=3;

double vectorsum(vector<double> arg, int size){
  double ans = 0;
  for(int i=0;i<size;i++){
    ans += (double)arg[i];
  }
  return ans;
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

void print_pr(){
    print("========================");
    for(int i=0;i<=max_node_num;i++){
        print(i<<": "<<pr[i]);
    }

}

string check_pr_status(vector<double> old_pr){
    double gsum = 0;
    for(int i=0;i<max_node_num+1;i++){
      auto tt = graph.find(i);
      int node = tt->first;
      auto outlinks = tt->second;
      double sum = 0;
      for(int i=0;i<outlinks.size();i++){
          if(outlinks[i]==-1) continue;
          auto it = graph.find(outlinks[i]);
          assert(it!=graph.end());
          sum += old_pr[outlinks[i]]/(it->second.size());
      }
      sum = alpha*sum + (1-alpha)/(max_node_num+1);
      gsum += sum;
      print(node<<" - Calculated: "<<sum<<" - In table: "<<pr[node]);
    }
    print("gsum="<<gsum);
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
            if(l.size()==1 && l[0]==-1){
              runtime.emit_intermediate(-2, pr[i]);
            }
            else
            {
              for(int j=0;j<l.size();j++){

                  if(l[j]==-1){break;}

                  assert(l[j]<=max_node_num);

                  runtime.emit_intermediate(l[j], pr[i]/l.size());
              }
            }
            runtime.emit_intermediate(i,0);
        }
    }
};

struct reduce_task : public mapreduce::reduce_task<int, double>
{
    template<typename Runtime, typename It>
    void operator()(Runtime &runtime, key_type const &key, It it, It ite) const
    {
        // print("YOLO");
        value_type sum = 0;
        for(auto itr=it;itr!=ite;itr++){
            sum += (double)(*itr);
        }
        if(key==-2){
            de = alpha * sum/(double)(max_node_num+1);
        }
        else
        {
          assert(key<=max_node_num);

          runtime.emit(key, (double)((1-alpha)+ alpha*(sum)));
          // runtime.emit(key, ((1-alpha)/(double)(max_node_num+1)+ alpha*(sum)/(double)(max_node_num+1)));
        }
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
    // print("heyo");
    mapreduce::specification spec;

    string filename = "diamond.txt";
    if (argc > 1)
        filename = string(argv[1]);

    if(argv[2][0]!='-' or argv[2][1]!='o'){
      cerr<<"Wrong Format"<<endl;
      cerr<<argv[2]<<endl;
    }

    // Reading file
    // print("File reading starts: "<<filename);
    ifstream fin;
    fin.open(filename);
    max_node_num = 0;
    int maxs = 0;
    int size = 0;
    if(fin){
        int s;
        while(fin>>s){
            size += 1;
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

    // print("File ready");

    for(int i=0;i<=max_node_num;i++){
        if(graph.find(i)==graph.end()) graph[i].push_back(-1);
    }

    // print(max_node_num);

    // print("Performing pre-processing");
    for(int i=0;i<max_node_num+1;i++){

      // Printing the pre-processing status
      // if((int)((i+1)/(size+1)*100)%5==0){
      //   print((float)(i+1)/(max_node_num+1)*100<<"% pre - processing done.");
      //   if((float)(i+1)/(max_node_num+1)*100==50){
      //     print("((((( HALFWAY THERE )))))");
      //   }
      // }

      for(int j=0;j<max_node_num+1;j++){
        auto itr = graph.find(j);
        vector<int> outlinks_to_j = itr->second;
        if(outlinks_to_j.size()==1 && outlinks_to_j[0] == -1){
          mm[i][j] = (alpha) * 1/(double)(max_node_num);
        }
        else
        {
          bool flag = -1;
          for(int k=0;k<outlinks_to_j.size();k++){
            if(outlinks_to_j[k]==i){
              mm[i][j] = (alpha) * 1/outlinks_to_j.size();
              flag = 1;
              break;
            }
          }
          if(flag==-1) mm[i][j] = 0;
        }
        mm[i][j] += (1-alpha) * (1/(double)max_node_num);
      }
    }

    // Setting up dangling vector
    // print("Danglings started");
    for(int i=0;i<max_node_num+1;i++){
      auto itr = graph.find(i);
      if(itr==graph.end()){
        print("Error(Dangling pointers overflow)");
      }
      auto ll = itr->second;
      if(ll.size()==0 && ll[0]==-1){
        dangling[i]=1;
      }
      else
      {
        dangling[i] = 0;
      }
    }
    de = 0;

    // Initializing prob distribution
    for(int i=0;i<=max_node_num;i++){
        pr[i] = (double)1/(max_node_num+1);
    }
    // print("Probabilities set up");

    // Specifying number of map tasks and reduce tasks
    int gtask = 8;
    if(max_node_num<=5){
      gtask = 1;
    }
    spec.map_tasks = gtask;
    int reduce_tasks = gtask;
    spec.reduce_tasks = reduce_tasks;
    pagerank::job::datasource_type number_source(0, max_node_num, (max_node_num+1)/reduce_tasks);
    pagerank::job job(number_source, spec);
    mapreduce::results result;

    // main ops started
    int iter = 0;
    double sum_all_pr = 0;

    double old_pr[max_node_num+1];
    for(int i=0;i<=max_node_num;i++){
        old_pr[i] = (double)(1/(double)(max_node_num+1));
    }
    double suma = 1e5;

    auto start = high_resolution_clock::now();

    while((abs(suma)>ACC_TOLERANCE) && iter<MAX_ITERATIONS)
    // for(int rr=0;rr<100;rr++)
    {
        iter += 1;
        // print(iter);
        #ifdef _DEBUG
            job.run<mapreduce::schedule_policy::sequential<pagerank::job> >(result);
        #else
            job.run<mapreduce::schedule_policy::cpu_parallel<pagerank::job> >(result);
        #endif

        for (auto it=job.begin_results(); it!=job.end_results(); ++it){
            pr[it->first] = it->second;
        }

        for(int i=0;i<=max_node_num;i++){
          pr[i] += (double)(pr[i]+de);
        }

        // double last_sum = 0;
        // for(int i=0;i<=max_node_num;i++){
        //   last_sum += pr[i];
        // }

        for(int a = 0;a<max_node_num+1;a++){
          double suu = 0;
          for(int b=0;b<max_node_num+1;b++){
            suu+=mm[a][b]*pr[b];
          }
          pr[a] = suu;
        }

        // for(int i=0;i<=max_node_num;i++){
        //   pr[i] = pr[i]/vectorsum(pr, max_node_num+1);
        // }

        suma = 0;
        for(int k=0;k<=max_node_num;k++){
            suma += (old_pr[k]-pr[k]);
        }

        for(int k=0;k<=max_node_num;k++){
            old_pr[k] = pr[k];
        }
    }

    // print("Num iterations: "<<iter);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    double last_sum = vectorsum(pr, max_node_num+1);
    for(int i=0;i<=max_node_num;i++){
      pr[i] = pr[i]/(double)last_sum;
    }

    ofstream outfile;
    outfile.open(argv[3]);

    double finalsum = vectorsum(pr, max_node_num+1);
    for(int i=0;i<=max_node_num;i++){
        outfile<<i<<" = "<<(double)pr[i]<<endl;
        // cout<<i<<" = "<<(double)pr[i]<<endl;
    }
    outfile<<"s = "<<finalsum<<endl;
    // cout<<"s = "<<finalsum<<endl;
    outfile.close();
    // print("Pagerank written successfully");

    // print("** Time of execution = "<<duration.count()/(double)1e6<<" **");
    print(filename<<", "<<duration.count()/(double)1e6);
    return 0;
}
