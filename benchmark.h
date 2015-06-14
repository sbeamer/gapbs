// See LICENSE.txt for license details.

#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <algorithm>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include "builder.h"
#include "graph.h"
#include "print_util.h"
#include "reader.h"
#include "timer.h"
#include "writer.h"


typedef int NodeID;
typedef NodeID WeightT;
typedef NodeWeight<NodeID, WeightT> WNode;

typedef CSRGraph<NodeID> Graph;
typedef CSRGraph<NodeID, WNode> WGraph;

typedef BuilderBase<NodeID, NodeID, NodeID> Builder;
typedef BuilderBase<NodeID, WNode, WeightT> WeightedBuilder;

typedef WriterBase<NodeID, NodeID> Writer;
typedef WriterBase<NodeID, WNode> WeightedWriter;


template<typename GraphT_>
class SourcePicker {
 public:
  SourcePicker(const GraphT_ &g, NodeID given_source=-1)
      : given_source(given_source), rng(8), udist(0, g.num_nodes()-1), g_(g) {}

  NodeID PickNext() {
    if (given_source != -1)
      return given_source;
    NodeID source;
    do {
      source = udist(rng);
    } while (g_.out_degree(source) == 0);
    return source;
  }

 private:
  NodeID given_source;
  std::mt19937 rng;
  std::uniform_int_distribution<NodeID> udist;
  const GraphT_ &g_;
};


template<typename KeyT, typename ValT>
std::vector<std::pair<ValT, KeyT>> TopK(
    const std::vector<std::pair<KeyT, ValT>> &to_sort, size_t k) {
  std::vector<std::pair<ValT, KeyT>> top_k;
  ValT min_so_far = 0;
  for (auto kvp : to_sort) {
    if ((top_k.size() < k) || (kvp.second > min_so_far)) {
      top_k.push_back(std::make_pair(kvp.second, kvp.first));
      std::sort(top_k.begin(), top_k.end(),
                std::greater<std::pair<ValT, KeyT>>());
      if (top_k.size() > k)
        top_k.resize(k);
      min_so_far = top_k.back().first;
    }
  }
  return top_k;
}


template<typename GraphT_, typename GraphFunc, typename AnalysisFunc>
void BenchmarkFunc(CLApp &cli, GraphT_ &g, GraphFunc f, AnalysisFunc a) {
  g.PrintStats();
  double search_total = 0;
  Timer search_timer;
  for (int iter=0; iter < cli.num_trials(); iter++) {
    search_timer.Start();
    auto result = f(g);
    search_timer.Stop();
    PrintTime("Trial Time", search_timer.Seconds());
    if (cli.do_analysis() && (iter == (cli.num_trials()-1))) {
      a(g, result);
    }
    search_total += search_timer.Seconds();
  }
  PrintTime("Search Time", search_total / cli.num_trials());
}

#endif  // BENCHMARK_H_
