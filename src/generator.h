// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <limits>
#include <random>

#include "graph.h"
#include "pvector.h"
#include "util.h"


/*
GAP Benchmark Suite
Class:  Generator
Author: Scott Beamer

Given scale and degree, generates edgelist for synthetic graph
 - Intended to be called from Builder
 - GenerateEL(uniform) generates and returns the edgelist
 - Can generate uniform random (uniform=true) or R-MAT graph according
   to Graph500 parameters (uniform=false)
 - Can also randomize weights within a weighted edgelist (InsertWeights)
 - Blocking/reseeding is for parallelism with deterministic output edgelist
*/


// maps to range [0,max_value], tailored to STL-style RNG
template <typename NodeID_, typename rng_t_,
          typename uNodeID_ = typename std::make_unsigned<NodeID_>::type>
class UniDist {
 public:
  UniDist(NodeID_ max_value, rng_t_ &rng): rng_(rng) {
    no_mod_ = rng_.max() == static_cast<uNodeID_>(max_value);
    mod_ = max_value + 1;
    uNodeID_ remainder_sub_1 = rng_.max() % mod_;
    if (remainder_sub_1 == mod_ - 1)
      cutoff_ = 0;
    else
      cutoff_ = rng_.max() - remainder_sub_1;
  }

  NodeID_ operator()() {
    uNodeID_ rand_num = rng_();
    if (no_mod_)
      return rand_num;
    if (cutoff_ != 0) {
      while (rand_num >= cutoff_)
        rand_num = rng_();
    }
    return rand_num % mod_;
  }

 private:
  rng_t_ &rng_;
  bool no_mod_;
  uNodeID_ mod_;
  uNodeID_ cutoff_;
};


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_,
          typename uNodeID_ = typename std::make_unsigned<NodeID_>::type,
          int uNodeID_bits_ = std::numeric_limits<uNodeID_>::digits,
          typename rng_t_ = typename std::conditional<(uNodeID_bits_ == 32),
                                                      std::mt19937,
                                                      std::mt19937_64>::type >
class Generator {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> WEdge;
  typedef pvector<Edge> EdgeList;

 public:
  Generator(int scale, int degree) {
    scale_ = scale;
    num_nodes_ = 1l << scale;
    num_edges_ = num_nodes_ * degree;
    if (num_nodes_ > std::numeric_limits<NodeID_>::max()) {
      std::cout << "NodeID type (max: " << std::numeric_limits<NodeID_>::max();
      std::cout << ") too small to hold " << num_nodes_ << std::endl;
      std::cout << "Recommend changing NodeID (typedef'd in src/benchmark.h)";
      std::cout << " to a wider type and recompiling" << std::endl;
      std::exit(-31);
    }
  }

  void PermuteIDs(EdgeList &el) {
    pvector<NodeID_> permutation(num_nodes_);
    rng_t_ rng(kRandSeed);
    #pragma omp parallel for
    for (NodeID_ n=0; n < num_nodes_; n++)
      permutation[n] = n;
    shuffle(permutation.begin(), permutation.end(), rng);
    #pragma omp parallel for
    for (int64_t e=0; e < num_edges_; e++)
      el[e] = Edge(permutation[el[e].u], permutation[el[e].v]);
  }

  EdgeList MakeUniformEL() {
    EdgeList el(num_edges_);
    #pragma omp parallel
    {
      rng_t_ rng;
      UniDist<NodeID_, rng_t_> udist(num_nodes_-1, rng);
      #pragma omp for
      for (int64_t block=0; block < num_edges_; block+=block_size) {
        rng.seed(kRandSeed + block/block_size);
        for (int64_t e=block; e < std::min(block+block_size, num_edges_); e++) {
          el[e] = Edge(udist(), udist());
        }
      }
    }
    return el;
  }

  EdgeList MakeRMatEL() {
    const uint32_t max = std::numeric_limits<uint32_t>::max();
    const uint32_t A = 0.57*max, B = 0.19*max, C = 0.19*max;
    EdgeList el(num_edges_);
    #pragma omp parallel
    {
      std::mt19937 rng;
      #pragma omp for
      for (int64_t block=0; block < num_edges_; block+=block_size) {
        rng.seed(kRandSeed + block/block_size);
        for (int64_t e=block; e < std::min(block+block_size, num_edges_); e++) {
          NodeID_ src = 0, dst = 0;
          for (int depth=0; depth < scale_; depth++) {
            uint32_t rand_point = rng();
            src = src << 1;
            dst = dst << 1;
            if (rand_point < A+B) {
              if (rand_point > A)
                dst++;
            } else {
              src++;
              if (rand_point > A+B+C)
                dst++;
            }
          }
          el[e] = Edge(src, dst);
        }
      }
    }
    PermuteIDs(el);
    // TIME_PRINT("Shuffle", std::shuffle(el.begin(), el.end(),
    //                                    std::mt19937()));
    return el;
  }

  EdgeList GenerateEL(bool uniform) {
    EdgeList el;
    Timer t;
    t.Start();
    if (uniform)
      el = MakeUniformEL();
    else
      el = MakeRMatEL();
    t.Stop();
    PrintTime("Generate Time", t.Seconds());
    return el;
  }

  static void InsertWeights(pvector<EdgePair<NodeID_, NodeID_>> &el) {}

  // Overwrites existing weights with random from [1,255]
  static void InsertWeights(pvector<WEdge> &el) {
    #pragma omp parallel
    {
      rng_t_ rng;
      UniDist<WeightT_, rng_t_> udist(254, rng);
      int64_t el_size = el.size();
      #pragma omp for
      for (int64_t block=0; block < el_size; block+=block_size) {
        rng.seed(kRandSeed + block/block_size);
        for (int64_t e=block; e < std::min(block+block_size, el_size); e++) {
          el[e].v.w = static_cast<WeightT_>(udist()+1);
        }
      }
    }
  }

 private:
  int scale_;
  int64_t num_nodes_;
  int64_t num_edges_;
  static const int64_t block_size = 1<<18;
};

#endif  // GENERATOR_H_
