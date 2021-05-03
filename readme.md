# Well Separated Pairs and the Traveling Salesman Problem

### Instructions for Running
>python tsp-nnp.py \<points file\> <separation factor> <flags:{-d, -bf}> \
>python tsp-nnp.py att48.tsp 5
>python tsp-nnp.py vlsi/xqf131.tsp 5

```
flags:
-d: debug info quadtree and WSP
-bf: brute force, turns off WSP for testing complexity without optimization
```

### Directory Structure
`tsp-bfp.py` - Brute Force + WSP Pruning algorithm \
`tsp-nnp.py` - Nearest Neighbor + WSP Pruning algorithm \
Point set data is located inside the `/data` directory \
WSP implementation inside `/wsp` directory

### Point data
We use the .TSP format common among TSP datasets.
