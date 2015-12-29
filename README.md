GAP Benchmark Suite [![Build Status](https://travis-ci.org/sbeamer/gapbs.svg)](https://travis-ci.org/sbeamer/gapbs)
===================

This is a pre-release version of the reference code for the upcoming [GAP](http://gap.cs.berkeley.edu/) [Benchmark Suite](http://gap.cs.berkeley.edu/benchmark.html). It is designed to be a portable high-performance baseline that only requires a compiler with support for C++11. For parallelism it uses OpenMP, but it can also be compiled to run serially without OpenMP. The details of the benchmark can be found in the [draft specification](http://arxiv.org/abs/1508.03619).

Kernels Included
----------------
+ Breadth-First Search (BFS)
+ Single-Source Shortest Paths (SSSP)
+ PageRank (PR)
+ Connected Components (CC)
+ Betweenness Centrality (BC)
+ Triangle Counting (TC)


Quick Start
-----------

Build the project:

    $ make

Override the default C++ compiler:

    $ CXX=g++-4.9 make

Test the build:

    $ make test

Run BFS on 1,024 vertices for 1 iteration:

    $ ./bfs -g 10 -n 1

Additional command line flags can be found with `-h`


Graph Loading
-------------

All of the binaries use the same command-line options for loading graphs:
+ `-g 20` generates a Kronecker graph with 2^20 vertices (Graph500 specifications)
+ `-u 20` generates a uniform random graph with 2^20 vertices (degree 16)
+ `-f graph.el` loads graph from file graph.el
+ `-sf graph.el` symmetrizes graph loaded from file graph.el

The graph loading infrastructure understands the following formats:
+ `.el` plain-text edge-list with an edge per line as _node1_ _node2_
+ `.wel` plain-text weighted edge-list with an edge per line as _node1_ _node2_ _weight_
+ `.gr` [9th DIMACS Implementation Challenge](http://www.dis.uniroma1.it/challenge9/download.shtml) format
+ `.graph` Metis format (used in [10th DIMACS Implementation Challenge](http://www.cc.gatech.edu/dimacs10/index.shtml))
+ `.mtx` [Matrix Market](http://math.nist.gov/MatrixMarket/formats.html) format
+ `.sg` serialized pre-built graph (use `converter` to make)
+ `.wsg` weighted serialized pre-built graph (use `converter` to make)


Future Features
---------------

+ Release integration support for OpenTuner
+ Scripts to perform official benchmark runs



How to Cite
-----------

Please cite this code by the benchmark specification:

Scott Beamer, Krste Asanović, David Patterson. [*The GAP Benchmark Suite*](http://arxiv.org/abs/1508.03619). arXiv:1508.03619 [cs.DC], 2015.
