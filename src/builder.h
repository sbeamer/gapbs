// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BUILDER_H_
#define BUILDER_H_

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <functional>
#include <type_traits>
#include <utility>

#include "command_line.h"
#include "generator.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "reader.h"
#include "timer.h"
#include "util.h"


/*
GAP Benchmark Suite
Class:  BuilderBase
Author: Scott Beamer

Given arguments from the command line (cli), returns a built graph
 - MakeGraph() will parse cli and obtain edgelist to call
   MakeGraphFromEL(edgelist) to perform the actual graph construction
 - edgelist can be from file (Reader) or synthetically generated (Generator)
 - Common case: BuilderBase typedef'd (w/ params) to be Builder (benchmark.h)
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_, bool invert = true>
class BuilderBase {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;

  const CLBase &cli_;
  bool symmetrize_;
  bool needs_weights_;
  bool in_place_ = false;
  int64_t num_nodes_ = -1;

 public:
  explicit BuilderBase(const CLBase &cli) : cli_(cli) {
    symmetrize_ = cli_.symmetrize();
    needs_weights_ = !std::is_same<NodeID_, DestID_>::value;
    in_place_ = cli_.in_place();
    if (in_place_ && needs_weights_) {
      std::cout << "In-place building (-m) does not support weighted graphs"
                << std::endl;
      exit(-30);
    }
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeID_> e) {
    return e.u;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> e) {
    return NodeWeight<NodeID_, WeightT_>(e.u, e.v.w);
  }

  NodeID_ FindMaxNodeID(const EdgeList &el) {
    NodeID_ max_seen = 0;
    #pragma omp parallel for reduction(max : max_seen)
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      max_seen = std::max(max_seen, e.u);
      max_seen = std::max(max_seen, (NodeID_) e.v);
    }
    return max_seen;
  }

  pvector<NodeID_> CountDegrees(const EdgeList &el, bool transpose) {
    pvector<NodeID_> degrees(num_nodes_, 0);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        fetch_and_add(degrees[e.u], 1);
      if ((symmetrize_ && !in_place_) || (!symmetrize_ && transpose))
        fetch_and_add(degrees[(NodeID_) e.v], 1);
    }
    return degrees;
  }

  static
  pvector<SGOffset> PrefixSum(const pvector<NodeID_> &degrees) {
    pvector<SGOffset> sums(degrees.size() + 1);
    SGOffset total = 0;
    for (size_t n=0; n < degrees.size(); n++) {
      sums[n] = total;
      total += degrees[n];
    }
    sums[degrees.size()] = total;
    return sums;
  }

  static
  pvector<SGOffset> ParallelPrefixSum(const pvector<NodeID_> &degrees) {
    const size_t block_size = 1<<20;
    const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
    pvector<SGOffset> local_sums(num_blocks);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset lsum = 0;
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++)
        lsum += degrees[i];
      local_sums[block] = lsum;
    }
    pvector<SGOffset> bulk_prefix(num_blocks+1);
    SGOffset total = 0;
    for (size_t block=0; block < num_blocks; block++) {
      bulk_prefix[block] = total;
      total += local_sums[block];
    }
    bulk_prefix[num_blocks] = total;
    pvector<SGOffset> prefix(degrees.size() + 1);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset local_total = bulk_prefix[block];
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++) {
        prefix[i] = local_total;
        local_total += degrees[i];
      }
    }
    prefix[degrees.size()] = bulk_prefix[num_blocks];
    return prefix;
  }

  // Removes self-loops and redundant edges
  // Side effect: neighbor IDs will be sorted
  void SquishCSR(const CSRGraph<NodeID_, DestID_, invert> &g, bool transpose,
                 DestID_*** sq_index, DestID_** sq_neighs) {
    pvector<NodeID_> diffs(g.num_nodes());
    DestID_ *n_start, *n_end;
    #pragma omp parallel for private(n_start, n_end)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose) {
        n_start = g.in_neigh(n).begin();
        n_end = g.in_neigh(n).end();
      } else {
        n_start = g.out_neigh(n).begin();
        n_end = g.out_neigh(n).end();
      }
      std::sort(n_start, n_end);
      DestID_ *new_end = std::unique(n_start, n_end);
      new_end = std::remove(n_start, new_end, n);
      diffs[n] = new_end - n_start;
    }
    pvector<SGOffset> sq_offsets = ParallelPrefixSum(diffs);
    *sq_neighs = new DestID_[sq_offsets[g.num_nodes()]];
    *sq_index = CSRGraph<NodeID_, DestID_>::GenIndex(sq_offsets, *sq_neighs);
    #pragma omp parallel for private(n_start)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose)
        n_start = g.in_neigh(n).begin();
      else
        n_start = g.out_neigh(n).begin();
      std::copy(n_start, n_start+diffs[n], (*sq_index)[n]);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> SquishGraph(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    DestID_ **out_index, *out_neighs, **in_index, *in_neighs;
    SquishCSR(g, false, &out_index, &out_neighs);
    if (g.directed()) {
      if (invert)
        SquishCSR(g, true, &in_index, &in_neighs);
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, in_index,
                                                in_neighs);
    } else {
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs);
    }
  }

  /*
  In-Place Graph Building Steps
    - sort edges and squish (remove self loops and redundant edges)
    - overwrite EdgeList's memory with outgoing neighbors
    - if graph not being symmetrized
      - finalize structures and make incoming structures if requested
    - if being symmetrized
      - search for needed inverses, make room for them, add them in place
  */
  void MakeCSRInPlace(EdgeList &el, DestID_*** index, DestID_** neighs,
                      DestID_*** inv_index, DestID_** inv_neighs) {
    // preprocess EdgeList - sort & squish in place
    std::sort(el.begin(), el.end());
    auto new_end = std::unique(el.begin(), el.end());
    el.resize(new_end - el.begin());
    auto self_loop = [](Edge e){ return e.u == e.v; };
    new_end = std::remove_if(el.begin(), el.end(), self_loop);
    el.resize(new_end - el.begin());
    // analyze EdgeList and repurpose it for outgoing edges
    pvector<NodeID_> degrees = CountDegrees(el, false);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    pvector<NodeID_> indegrees = CountDegrees(el, true);
    *neighs = reinterpret_cast<DestID_*>(el.data());
    for (Edge e : el)
      (*neighs)[offsets[e.u]++] = e.v;
    size_t num_edges = el.size();
    el.leak();
    // revert offsets by shifting them down
    for (NodeID_ n = num_nodes_; n >= 0; n--)
      offsets[n] = n != 0 ? offsets[n-1] : 0;
    if (!symmetrize_) {   // not going to symmetrize so no need to add edges
      size_t new_size = num_edges * sizeof(DestID_);
      *neighs = static_cast<DestID_*>(std::realloc(*neighs, new_size));
      *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs);
      if (invert) {       // create inv_neighs & inv_index for incoming edges
        pvector<SGOffset> inoffsets = ParallelPrefixSum(indegrees);
        *inv_neighs = new DestID_[inoffsets[num_nodes_]];
        *inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inoffsets,
                                                          *inv_neighs);
        for (NodeID_ u = 0; u < num_nodes_; u++) {
          for (DestID_* it = (*index)[u]; it < (*index)[u+1]; it++) {
            NodeID_ v = static_cast<NodeID_>(*it);
            (*inv_neighs)[inoffsets[v]] = u;
            inoffsets[v]++;
          }
        }
      }
    } else {              // symmetrize graph by adding missing inverse edges
      // Step 1 - count number of needed inverses
      pvector<NodeID_> invs_needed(num_nodes_, 0);
      for (NodeID_ u = 0; u < num_nodes_; u++) {
        for (SGOffset i = offsets[u]; i < offsets[u+1]; i++) {
          DestID_ v = (*neighs)[i];
          bool inv_found = std::binary_search(*neighs + offsets[v],
                                              *neighs + offsets[v+1],
                                              static_cast<DestID_>(u));
          if (!inv_found)
            invs_needed[v]++;
        }
      }
      // increase offsets to account for missing inverses, realloc neighs
      SGOffset total_missing_inv = 0;
      for (NodeID_ n = 0; n < num_nodes_; n++) {
        offsets[n] += total_missing_inv;
        total_missing_inv += invs_needed[n];
      }
      offsets[num_nodes_] += total_missing_inv;
      size_t newsize = (offsets[num_nodes_] * sizeof(DestID_));
      *neighs = static_cast<DestID_*>(std::realloc(*neighs, newsize));
      if (*neighs == nullptr) {
        std::cout << "Call to realloc() failed" << std::endl;
        exit(-33);
      }
      // Step 2 - spread out existing neighs to make room for inverses
      //   copies backwards (overwrites) and inserts free space at starts
      SGOffset tail_index = offsets[num_nodes_] - 1;
      for (NodeID_ n = num_nodes_ - 1; n >= 0; n--) {
        SGOffset new_start = offsets[n] + invs_needed[n];
        for (SGOffset i = offsets[n+1]-1; i >= new_start; i--) {
          (*neighs)[tail_index] = (*neighs)[i - total_missing_inv];
          tail_index--;
        }
        total_missing_inv -= invs_needed[n];
        tail_index -= invs_needed[n];
      }
      // Step 3 - add missing inverse edges into free spaces from Step 2
      for (NodeID_ u = 0; u < num_nodes_; u++) {
        for (SGOffset i = offsets[u] + invs_needed[u]; i < offsets[u+1]; i++) {
          DestID_ v = (*neighs)[i];
          bool inv_found = std::binary_search(
                             *neighs + offsets[v] + invs_needed[v],
                             *neighs + offsets[v+1],
                             static_cast<DestID_>(u));
          if (!inv_found) {
            (*neighs)[offsets[v] + invs_needed[v] -1] = static_cast<DestID_>(u);
            invs_needed[v]--;
          }
        }
      }
      for (NodeID_ n = 0; n < num_nodes_; n++)
        std::sort(*neighs + offsets[n], *neighs + offsets[n+1]);
      *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs);
    }
  }

  /*
  Graph Building Steps (for CSR):
    - Read edgelist once to determine vertex degrees (CountDegrees)
    - Determine vertex offsets by a prefix sum (ParallelPrefixSum)
    - Allocate storage and set points according to offsets (GenIndex)
    - Copy edges into storage
  */
  void MakeCSR(const EdgeList &el, bool transpose, DestID_*** index,
               DestID_** neighs) {
    pvector<NodeID_> degrees = CountDegrees(el, transpose);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    *neighs = new DestID_[offsets[num_nodes_]];
    *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        (*neighs)[fetch_and_add(offsets[e.u], 1)] = e.v;
      if (symmetrize_ || (!symmetrize_ && transpose))
        (*neighs)[fetch_and_add(offsets[static_cast<NodeID_>(e.v)], 1)] =
            GetSource(e);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraphFromEL(EdgeList &el) {
    DestID_ **index = nullptr, **inv_index = nullptr;
    DestID_ *neighs = nullptr, *inv_neighs = nullptr;
    Timer t;
    t.Start();
    if (num_nodes_ == -1)
      num_nodes_ = FindMaxNodeID(el)+1;
    if (needs_weights_)
      Generator<NodeID_, DestID_, WeightT_>::InsertWeights(el);
    if (in_place_) {
      MakeCSRInPlace(el, &index, &neighs, &inv_index, &inv_neighs);
    } else {
      MakeCSR(el, false, &index, &neighs);
      if (!symmetrize_ && invert) {
        MakeCSR(el, true, &inv_index, &inv_neighs);
      }
    }
    t.Stop();
    PrintTime("Build Time", t.Seconds());
    if (symmetrize_)
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs);
    else
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs,
                                                inv_index, inv_neighs);
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraph() {
    CSRGraph<NodeID_, DestID_, invert> g;
    {  // extra scope to trigger earlier deletion of el (save memory)
      EdgeList el;
      if (cli_.filename() != "") {
        Reader<NodeID_, DestID_, WeightT_, invert> r(cli_.filename());
        if ((r.GetSuffix() == ".sg") || (r.GetSuffix() == ".wsg")) {
          return r.ReadSerializedGraph();
        } else {
          el = r.ReadFile(needs_weights_);
        }
      } else if (cli_.scale() != -1) {
        Generator<NodeID_, DestID_> gen(cli_.scale(), cli_.degree());
        el = gen.GenerateEL(cli_.uniform());
      }
      g = MakeGraphFromEL(el);
    }
    if (in_place_)
      return g;
    else
      return SquishGraph(g);
  }

  // Relabels (and rebuilds) graph by order of decreasing degree
  static
  CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    if (g.directed()) {
      std::cout << "Cannot relabel directed graph" << std::endl;
      std::exit(-11);
    }
    Timer t;
    t.Start();
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++)
      degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
    std::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
              std::greater<degree_node_p>());
    pvector<NodeID_> degrees(g.num_nodes());
    pvector<NodeID_> new_ids(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      degrees[n] = degree_id_pairs[n].first;
      new_ids[degree_id_pairs[n].second] = n;
    }
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
    #pragma omp parallel for
    for (NodeID_ u=0; u < g.num_nodes(); u++) {
      for (NodeID_ v : g.out_neigh(u))
        neighs[offsets[new_ids[u]]++] = new_ids[v];
      std::sort(index[new_ids[u]], index[new_ids[u]+1]);
    }
    t.Stop();
    PrintTime("Relabel", t.Seconds());
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
  }
};

#endif  // BUILDER_H_
