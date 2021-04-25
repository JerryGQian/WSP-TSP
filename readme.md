# Well Separated Pairs and the Traveling Salesman Problem

### Instructions for Running
>python tsp-bfp.py \<points file\> <separation factor> <flags:{-d, -bf}> \
>python tsp-bfp.py points1.txt 1 -d 

```
flags:
-d: debug info quadtree and WSP
-bf: brute force, turns off WSP for testing complexity without optimization
```

### Directory Structure
`tsp-bfp.py` - Brute Force + WSP Pruning algorithm \
Point set data is located inside the `/data` directory \
WSP implementation inside `/wsp` directory

### Point data
You can specify your own point data by creating a txt file in the `/data` directory. \
File format: \
`bounds:<xMin>,<yMin>,<xMax>,<yMax>` - bounds set the min and max dimensions of the point space. Ex: `bounds:0,0,100,100` \
`<int>,<int>` - point, should be within bounds. Ex: `25,30` \
`# comment` \
` ` - empty lines for organization acceptable
