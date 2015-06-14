// See LICENSE.txt for license details.

#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"

using namespace std;


typedef float ScoreT;
const float kDamp = 0.85;

pvector<ScoreT> PageRankPull(Graph &g, int num_iterations) {
  const ScoreT init_score = 1.0f / g.num_nodes();
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
  pvector<ScoreT> scores(g.num_nodes(), init_score);
  pvector<ScoreT> contrib(g.num_nodes(), 0);
  for (int iter=0; iter < num_iterations; iter++) {
    ScoreT error = 0;
    #pragma omp parallel for
    for (NodeID n=0; n < g.num_nodes(); n++)
      contrib[n] = scores[n] / g.out_degree(n);
    #pragma omp parallel for reduction(+ : error)
    for (NodeID u=0; u < g.num_nodes(); u++) {
      ScoreT sum = 0;
      for (NodeID v : g.in_neigh(u))
        sum += contrib[v];
      ScoreT old_score = scores[u];
      scores[u] = base_score + kDamp * sum;
      error += fabs(scores[u] - old_score);
    }
    cout << " " << iter << "    " << error << endl;
  }
  return scores;
}

void PrintTopScores(Graph &g, pvector<ScoreT> &scores) {
  vector<pair<NodeID, ScoreT>> score_pairs(g.num_nodes());
  for (NodeID n=0; n < g.num_nodes(); n++) {
    score_pairs[n] = make_pair(n, scores[n]);
  }
  int k = 5;
  vector<pair<ScoreT, NodeID>> top_k = TopK(score_pairs, k);
  k = min(k, static_cast<int>(top_k.size()));
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
}

int main(int argc, char* argv[]) {
  CLIterApp cli(argc, argv, "pagerank", 15);
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  auto PRBound = [&cli] (Graph &g) { return PageRankPull(g, cli.num_iters()); };
  BenchmarkFunc(cli, g, PRBound, PrintTopScores);
  return 0;
}
