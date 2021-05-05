# Well Separated Pairs and the Traveling Salesman Problem

### Instructions for Running
>python tsp-nnp.py \<TSP dataset file\> <separation factor> <flags:{-d, -bf}> \
>python tsp-nnp.py att48.tsp 5
>python tsp-nnp.py vlsi/xqf131.tsp 5
>python tsp-wgmst.py vlsi/xqf131.tsp 1 -pr

```
flags:
-d: debug info quadtree and WSP
-bf: brute force, turns off WSP for testing complexity without optimization
```

### Directory Structure
`tsp-bfp.py` - Brute Force + WSP Pruning algorithm - Exponential \
`tsp-nnp.py` - Nearest Neighbor + WSP Pruning algorithm - Polynomial \
`tsp-sps.py` - Subproblem Sorting with WSP - Polynomial - questionable results \
`tsp-wgmst.py` - WSP Graph MST, builds MST out of graph found from WSP min projections - In progress. Motivation: https://www.cs.umd.edu/~mount/Indep/Chaojun_Li/final-rept.pdf \
Point set data is located inside the `/data` directory \
WSP implementation inside `/wsp` directory

### Point data
We use the .TSP format common among TSP datasets.
