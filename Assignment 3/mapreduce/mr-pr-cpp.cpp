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
// #include <iostream>
#include <bits/stdc++.h>
#include "table.h"

#define print(x) cout<<x<<'\n'

namespace pagerank{

// std::pair<std::int, std::vector<std::int, std::int> > graph;
vector<pair<int, int> > inGraph;
vector<pair<int, vector<int> >  > outGraph = {{1, {2}}, {2, {1}}};

struct map_task : public mapreduce::map_task<int, vector<int> >
{
	template<typename Runtime>
	void operator()(Runtime &runtime, int const &key, vector<int> const &value) const
	{

	}
};

struct reduce_task : public mapreduce::reduce_task<bool, long>
{
    template<typename Runtime, typename It>
    void operator()(Runtime &runtime, int const &key, It it, It ite) const
    {

    }
};

typedef mapreduce::job<pagerank::map_task, pagerank::reduce_task, mapreduce::null_combiner> job;

}

int main(int argc, char const *argv[])
{
	print("Hello world");
	Table t;
	
	t.read_file("diamond.txt");
	t.print_table();
	
	cerr << "Calculating pagerank..." << endl;
    t.pagerank();
    cerr << "Done calculating!" << endl;
    
    return 0;
}