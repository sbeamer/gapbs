// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <functional>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "bitmap.h"
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


void PBFS(const Graph &g, NodeID source, pvector<NodeID> &path_counts,
    Bitmap &succ, vector<SlidingQueue<NodeID>::iterator> &depth_index,
    SlidingQueue<NodeID> &queue) {
  pvector<NodeID> depths(g.num_nodes(), -1);
  depths[source] = 0;
  path_counts[source] = 1;
  queue.push_back(source);
  depth_index.push_back(queue.begin());
  queue.slide_window();
  const NodeID* g_out_start = g.out_neigh(0).begin();
  #pragma omp parallel
  {
    int depth = 0;
    QueueBuffer<NodeID> lqueue(queue);
    while (!queue.empty()) {
      #pragma omp single
      depth_index.push_back(queue.begin());
      depth++;
      #pragma omp for
      for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
        NodeID u = *q_iter;
        for (NodeID &v : g.out_neigh(u)) {
          if ((depths[v] == -1) && (compare_and_swap(depths[v], -1, depth))) {
            lqueue.push_back(v);
          }
          if (depths[v] == depth) {
            succ.set_bit_atomic(&v - g_out_start);
            fetch_and_add(path_counts[v], path_counts[u]);
          }
        }
      }
      lqueue.flush();
      #pragma omp barrier
      #pragma omp single
      queue.slide_window();
    }
  }
  depth_index.push_back(queue.begin());
}


pvector<ScoreT> Brandes(const Graph &g, SourcePicker<Graph> &sp,
                        NodeID num_iters) {
  Timer t;
  t.Start();
  pvector<ScoreT> scores(g.num_nodes(), 0);
  pvector<NodeID> path_counts(g.num_nodes());
  Bitmap succ(g.num_edges_directed());
  vector<SlidingQueue<NodeID>::iterator> depth_index;
  SlidingQueue<NodeID> queue(g.num_nodes());
  t.Stop();
  PrintStep("a", t.Seconds());
  const NodeID* g_in_start = g.in_neigh(0).begin();
  for (NodeID iter=0; iter < num_iters; iter++) {
    NodeID source = sp.PickNext();
    cout << "source: " << source << endl;
    t.Start();
    path_counts.fill(0);
    depth_index.resize(0);
    queue.reset();
    succ.reset();
    PBFS(g, source, path_counts, succ, depth_index, queue);
    t.Stop();
    PrintStep("b", t.Seconds());
    pvector<ScoreT> deltas(g.num_nodes(), 0);
    t.Start();
    for (int d=depth_index.size()-2; d >= 0; d--) {
      #pragma omp parallel for
      for (auto it = depth_index[d]; it < depth_index[d+1]; it++) {
        NodeID u = *it;
        ScoreT delta_u = 0;
        for (NodeID &v : g.in_neigh(u)) {
          if (succ.get_bit(&v - g_in_start)) {
            delta_u += static_cast<ScoreT>(path_counts[u]) /
                       static_cast<ScoreT>(path_counts[v]) * (1 + deltas[v]);
          }
        }
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
  if (cli.num_iters() > 1 && cli.start_vertex() != -1) {
    cout << "Warning: iterating from same source (-r & -k)" << endl;
  }
  Builder b(cli);
  Graph g = b.MakeGraph();
  SourcePicker<Graph> sp(g, cli.start_vertex());
  auto BCBound =
    [&sp, &cli] (const Graph &g) { return Brandes(g, sp, cli.num_iters()); };
  BenchmarkKernel(cli, g, BCBound, PrintTopScores);
  return 0;
}
