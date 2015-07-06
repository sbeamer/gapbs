// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <functional>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "print_util.h"
#include "pvector.h"
#include "sliding_queue.h"
#include "timer.h"


using namespace std;
typedef float ScoreT;


void BFS(Graph &g, NodeID source, pvector<NodeID> &path_counts,
    pvector<NodeID> &succ, pvector<SGOffset> &succ_tails,
    vector<SlidingQueue<NodeID>::iterator> &depth_index,
    SlidingQueue<NodeID> &queue) {
  pvector<NodeID> depths(g.num_nodes(), -1);
  depths[source] = 0;
  path_counts[source] = 1;
  queue.Push(source);
  depth_index.push_back(queue.begin());
  queue.SlideWindow();
  int depth = 0;
  while (!queue.Empty()) {
    depth_index.push_back(queue.begin());
    depth++;
    #pragma omp parallel for
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
      NodeID u = *q_iter;
      for (NodeID v : g.out_neigh(u)) {
        if ((depths[v] == -1) && (compare_and_swap(depths[v], -1, depth))) {
          queue.Push(v);
        }
        if (depths[v] == depth) {
          succ[fetch_and_add(succ_tails[u], 1)] = v;
          path_counts[v] += path_counts[u];
        }
      }
    }
    queue.SlideWindow();
  }
  depth_index.push_back(queue.begin());
}


pvector<ScoreT> Brandes(Graph &g, NodeID source, NodeID num_iters) {
  cout << "source: " << source << endl;
  Timer t;
  t.Start();
  pvector<ScoreT> scores(g.num_nodes(), 0);
  pvector<NodeID> path_counts(g.num_nodes());
  pvector<NodeID> succ(g.num_edges_directed());
  pvector<SGOffset> succ_heads = g.VertexOffsets();
  pvector<SGOffset> succ_tails(succ_heads.begin(), succ_heads.end());
  vector<SlidingQueue<NodeID>::iterator> depth_index;
  SlidingQueue<NodeID> queue(g.num_nodes());
  t.Stop();
  PrintStep("a", t.Seconds());
  for (NodeID iter=0; iter < num_iters; iter++) {
    path_counts.fill(0);
    depth_index.resize(0);
    queue.Reset();
    t.Start();
    BFS(g, source, path_counts, succ, succ_tails, depth_index, queue);
    t.Stop();
    PrintStep("b", t.Seconds());
    pvector<ScoreT> deltas(g.num_nodes(), 0);
    t.Start();
    for (int d=depth_index.size()-2; d >= 0; d--) {
      #pragma omp parallel for
      for (auto it = depth_index[d]; it < depth_index[d+1]; it++) {
        NodeID u = *it;
        ScoreT delta_u = 0;
        for (SGOffset i=succ_heads[u]; i < succ_tails[u]; i++) {
          NodeID v = succ[i];
          delta_u += static_cast<ScoreT>(path_counts[u]) /
                     static_cast<ScoreT>(path_counts[v]) * (1 + deltas[v]);
        }
        succ_tails[u] = succ_heads[u];    // resets succ_tails
        deltas[u] = delta_u;
        scores[u] += delta_u;
      }
    }
    t.Stop();
    PrintStep("p", t.Seconds());
  }
  return scores;
}


void PrintTopScores(const Graph &g, const pvector<ScoreT> &scores) {
  vector<pair<NodeID, ScoreT>> score_pairs(g.num_nodes());
  for (NodeID n=0; n < g.num_nodes(); n++) {
    score_pairs[n] = make_pair(n, scores[n]);
  }
  int k = 5;
  vector<pair<ScoreT, NodeID>> top_k = TopK(score_pairs, k);
  if (!score_pairs.empty()) {
    ScoreT top_score = top_k[0].first;
    for (auto kvp : top_k)
      cout << kvp.second << ":" << kvp.first / top_score << endl;
  }
}


int main(int argc, char* argv[]) {
  CLIterApp cli(argc, argv, "betweenness-centrality", 1);
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  SourcePicker<Graph> sp(g, cli.start_vertex());
  auto BCBound =
    [&sp, &cli] (Graph &g) { return Brandes(g, sp.PickNext(),
                                            cli.num_iters()); };
  BenchmarkFunc(cli, g, BCBound, PrintTopScores);
  return 0;
}
