// See LICENSE.txt for license details.

#include <iostream>
#include <vector>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "reader.h"
#include "sliding_queue.h"
#include "timer.h"

using namespace std;


/*

Using optimization of precomputing degrees in bulk in beginning and storing
them in parent array as negative numbers. Thus the encoding of parent is:
parent[x] < 0 implies it is -out_degree(x)
parent[x] >= 0 implies it is parent(x)

*/



long BUStep(Graph &g, pvector<NodeID> &parent,
                   Bitmap &front, Bitmap &next) {
  long awake_count = 0;
  next.reset();
  #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic,1024)
  for (NodeID u=0; u < g.num_nodes(); u++) {
    if (parent[u] < 0) {
      for (NodeID v : g.in_neigh(u)) {
        if (front.get_bit(v)) {
          parent[u] = v;
          awake_count++;
          next.set_bit(u);
          break;
        }
      }
    }
  }
  return awake_count;
}


long TDStep(Graph &g, pvector<NodeID> &parent, SlidingQueue<NodeID> &queue) {
  long scout_count = 0;
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue);
    #pragma omp for reduction(+ : scout_count)
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
      NodeID u = *q_iter;
      for (NodeID v : g.out_neigh(u)) {
        NodeID curr_val = parent[v];
        if (curr_val < 0) {
          if (compare_and_swap(parent[v], curr_val, u)) {
            lqueue.Push(v);
            scout_count += -curr_val;
          }
        }
      }
    }
    lqueue.Flush();
  }
  return scout_count;
}


void QueueToBitmap(SlidingQueue<NodeID> &queue, Bitmap &bm) {
  #pragma omp parallel for
  for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
    NodeID u = *q_iter;
    bm.set_bit_atomic(u);
  }
}

void BitmapToQueue(Graph &g, Bitmap &bm, SlidingQueue<NodeID> &queue) {
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue);
    #pragma omp for
    for (NodeID n=0; n < g.num_nodes(); n++)
      if (bm.get_bit(n))
        lqueue.Push(n);
    lqueue.Flush();
  }
  queue.SlideWindow();
}

pvector<NodeID> InitParent(Graph &g) {
  pvector<NodeID> parent(g.num_nodes());
  #pragma omp parallel for
  for (NodeID n=0; n < g.num_nodes(); n++)
    parent[n] = g.out_degree(n)!=0 ? -g.out_degree(n) : -1;
  return parent;
}


pvector<NodeID> DOBFS(Graph &g, NodeID source, int alpha=26, int beta=72) {
  cout << "source: " << source << endl;
  Timer t;
  t.Start();
  pvector<NodeID> parent = InitParent(g);
  t.Stop();
  PrintStep("i", t.Seconds());
  parent[source] = source;
  SlidingQueue<NodeID> queue(g.num_nodes());
  queue.Push(source);
  queue.SlideWindow();
  Bitmap curr(g.num_nodes());
  curr.reset();
  Bitmap front(g.num_nodes());
  front.reset();
  long edges_to_check = g.num_edges_directed();
  long scout_count = g.out_degree(source);
  while (!queue.Empty()) {
    if (scout_count > edges_to_check / alpha) {
      long awake_count;
      TIME_OP(t, QueueToBitmap(queue, front));
      PrintStep("e", t.Seconds());
      queue.SlideWindow();
      do {
        t.Start();
        awake_count = BUStep(g, parent, front, curr);
        front.swap(curr);
        t.Stop();
        PrintStep("bu", t.Seconds(), awake_count);
      } while (awake_count > g.num_nodes() / beta);
      TIME_OP(t, BitmapToQueue(g, front, queue));
      PrintStep("c", t.Seconds());
      scout_count = 1;
    } else {
      t.Start();
      edges_to_check -= scout_count;
      scout_count = TDStep(g, parent, queue);
      queue.SlideWindow();
      t.Stop();
      PrintStep("td", t.Seconds(), queue.size());
    }
  }
  return parent;
}


void PrintBFSStats(Graph &g, pvector<NodeID> &bfs_tree) {
  long tree_size = 0;
  long n_edges = 0;
  for (NodeID n=0; n < g.num_nodes(); n++) {
    if (bfs_tree[n] >= 0) {
      n_edges += g.out_degree(n);
      tree_size++;
    }
  }
  cout << "BFS Tree has " << tree_size << " nodes and ";
  cout << n_edges << " edges" << endl;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  SourcePicker<Graph> sp(g, cli.start_vertex());
  auto BFSBound = [&sp] (Graph &g) { return DOBFS(g, sp.PickNext()); };
  BenchmarkFunc(cli, g, BFSBound, PrintBFSStats);
  return 0;
}
