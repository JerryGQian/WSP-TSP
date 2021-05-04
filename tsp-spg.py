from wsp import wsp
from wsp import ds
from wsp import util
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import time

# SPG = sub problem graph

# run algorithm
# >> python tsp-spg.py <points file> <separation factor> <quadtree:{-pr, -point/-p}> <flags:{-d, -bf}>
# >> python tsp-spg.py att48.tsp 1 -p -d
# -d: debug info quadtree and WSP
# -bf: brute force (turns off WSP)

timeInit = time.perf_counter()

filename = "data/points.txt"
s = 1           # default separation factor
wsp_mode = True # uses WSPs
debug = False   # debug info for Quadtree and WSPs
quadtree = ds.PointQuadTree

if (len(sys.argv) >= 2):
  filename = "data/" + sys.argv[1]
if (len(sys.argv) >= 3):
  s = int(sys.argv[2])
# check flags
for arg in sys.argv:
    if arg == "-pr":
        quadtree = ds.PRQuadTree
    if arg == "-point" or arg == "-p":
        quadtree = ds.PointQuadTree
    if arg == "-d":
        debug = True
    if arg == "-bf":
        wsp_mode = False
# build WSP tree
wspTreeNode = wsp.runWSP(filename, s, debug, quadtree)

timeStart = time.perf_counter()

# calculate well separated dictionary
points = wspTreeNode.get_points()

points = wspTreeNode.get_points()
num_points = len(points)

graph = {}
for p in points:
    graph[p] = set()

def build_wsp_graph(cur_node):
    if len(cur_node.connection) > 0:
        for node in cur_node.connection:
            p1, p2 = ds.min_proj(cur_node, node)
            graph[p1].add(p2)
            graph[p2].add(p1)

    if cur_node.divided:
        build_wsp_graph(cur_node.ne)
        build_wsp_graph(cur_node.nw)
        build_wsp_graph(cur_node.sw)
        build_wsp_graph(cur_node.se)

wsp_count = 0
build_wsp_graph(wspTreeNode)
print("graph:", graph)

edge_count_graph = set()
for p1 in graph:
    for p2 in graph[p1]:
        if (p1, p2) not in edge_count_graph:
            edge_count_graph.add((p1, p2))
            edge_count_graph.add((p2, p1))
            wsp_count += 1
            #wsp.ax[1].plot([p1.x, p2.x],[p1.y, p2.y], color="red")

print("___________________________________________________________")
print(num_points, "points")
print(int(wsp_count/2), "WSPs found")
print(f"Loaded in {timeStart - timeInit:0.4f} seconds")


queue = [(points[0], 0, None)]
added = set()
mst = {}
#min_span_tree_cost = 0
# Prims alg
while queue:
    # Choose the adjacent node with the least edge cost
    minCost = float('inf')
    minEdge = None
    for tup in queue:
        if tup[1] < minCost:
            minCost = tup[1]
            minEdge = tup
    queue.remove(minEdge)
    node = minEdge[0]
    prev = minEdge[2]
    #print(node)
    if prev != None and node not in added:
        if prev not in mst:
            mst[prev] = set()
        
        mst[prev].add(node)
    #cost = queue[node]

    if node not in added:
        #min_span_tree_cost += cost
        added.add(node)
        #print("neighbor:", graph[node])
        for neighbor in graph[node]:
            if neighbor not in added:
                queue.append((neighbor, node.distance_to(neighbor), node))

print("mst:", mst)
drawn = set()
def draw_mst(node):
    #print(node)
    if node in mst and node not in drawn:
        for neighbor in mst[node]:
            wsp.ax[1].plot([node.x, neighbor.x],[node.y, neighbor.y], color="blue")
            if neighbor in mst:
                draw_mst(neighbor)
            drawn.add(node)
print(points[0])
draw_mst(points[0])

path = []
# traverse MST
def traverse_mst_for_path(prev, cur):
    def calc_angle(node):
        x = 0

    # order paths
    ordered_children = []
    if cur in mst:
        for child in mst[cur]:
            traverse_mst_for_path(cur, child)

    path.append(cur)

traverse_mst_for_path(None, points[0])

# find shortest permutation
minSolution = path
minDist = util.calcDist(path)

timeEnd = time.perf_counter()

for i in range(len(minSolution) - 1):
    wsp.ax[1].plot([minSolution[i].x, minSolution[i+1].x],[minSolution[i].y, minSolution[i+1].y], color="red")
wsp.ax[0].set_title(f"#WSP={wsp_count}, s={s}")
wsp.ax[1].set_title(f"TSP Path: n={len(points)}, length={minDist:0.4f}")

print("")
print("Solution:", minSolution)
print("Solution Distance:", minDist)
print(1, "permutations examined")
print(f"Solution found in {timeEnd - timeStart:0.4f} seconds")
print("___________________________________________________________")
plt.show()