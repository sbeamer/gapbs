GAP Benchmark Suite
===================

This is a pre-release copy of the reference code for the upcoming [GAP](https://www.cs.berkeley.edu/~sbeamer/gap/) Benchmark Suite .


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

Build the project with
    $ make

Run BFS on 1,000 vertices for 1 iteration
    $ ./bfs -g 10 -n 1

Additional command line flags can be found with `-h`


Graph Loading
-------------

All of the binaries use the same command-line options for loading graphs.
+ `-g 20` generates a SCALE=20 Kronecker graph (Graph500 specifications)
+ `-u 20` generates a SCALE=20 uniform random graph (degree 16)
+ `-f graph.el` loads graph from file graph.el
+ `-sf graph.el` symmetrizes graph loaded from file graph.el 

The graph loading infrastructure understands the following formats:
+ `.el` plain-text edge-list with an edge per line as <node1> <node2>
+ `.wel` plain-text weighted edge-list with an edge per line as <node1> <node2> <weight>
+ `.gr` DIMACS Grand Challenge format
+ `.sg` serialized pre-built graph (use converter to make)
+ `.sg` weighted serialized pre-built graph (use converter to make)
