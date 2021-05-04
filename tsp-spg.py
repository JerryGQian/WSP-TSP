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
print(graph)

edge_count_graph = set()
for p1 in graph:
    for p2 in graph[p1]:
        if (p1, p2) not in edge_count_graph:
            edge_count_graph.add((p1, p2))
            edge_count_graph.add((p2, p1))
            wsp_count += 1
            wsp.ax[1].plot([p1.x, p2.x],[p1.y, p2.y], color="red")

print("___________________________________________________________")
print(num_points, "points")
print(int(wsp_count/2), "WSPs found")
print(f"Loaded in {timeStart - timeInit:0.4f} seconds")

'''
# traversal
def closest_sub(p1, sublist):
  """Min dist point of sublist from p1"""
  mind = 99999999
  minPoint = None
  for p2 in sublist:
    dist = p1.distance_to(p2)
    if dist < mind:
        mind = dist
        minPoint = p2
  return minPoint, mind

def find_subpath(src, glist):
    if src == None: # try all points
        minPath = []
        minLen = float('inf')
        for p in glist:
            rem = glist.copy()
            rem.remove(p)
            #print("new subpath", p)
            if isinstance(p, list):
                path = find_subpath(None, p)
            else:
                path = find_path(p, rem)
            length = util.calcDist(path)
            if length < minLen:
                minLen = length
                minPath = path
        return path
    else:
        start,_ = closest_sub(src, glist)
        rem = glist.copy()
        rem.remove(start)
        #print("srced subpath", start)
        return find_path(start, rem)

def find_subpath_ends(src, dest, glist):
    if src == None: # try all points
        minPath = []
        minLen = float('inf')
        for p in glist:
            rem = glist.copy()
            rem.remove(p)
            #print("new subpath", p)
            if isinstance(p, list):
                path = find_subpath(None, p)
            else:
                path = find_path(p, rem)
            length = util.calcDist(path)
            if length < minLen:
                minLen = length
                minPath = path
        return path
    else:
        start,_ = closest_sub(src, glist)
        rem = glist.copy()
        rem.remove(start)
        #print("srced subpath", start)
        return find_path(start, rem)

def find_path(start, rem):
    path = [start]
    if len(path) == num_points:
        return [path]
    while len(rem) > 0:
        last_point = path[len(path) - 1]
        #print(last_point)
        minNext = None
        minNextDist = float('inf')
        for r in rem: # look through rem for next
            if isinstance(r, list):
                _, curNextDist = closest_sub(last_point, r)
                #path += find_subpath(last_point, r)
                #rem.remove(r)
            else:
                #print(last_point)
                curNextDist = last_point.distance_to(r)
            if curNextDist < minNextDist:
                minNext = r
                minNextDist = curNextDist
        if minNext == None:
            minNext = rem[0] # shouldnt be needed?
        if isinstance(minNext, list):
            path += find_subpath(last_point, minNext)
            #print(minNext)
        else:
            #print("should be point",minNext, len(rem))
            path.append(minNext)
        rem.remove(minNext)
    return path

# search for permutations
perms = []
for item in grouped_points:
    if isinstance(item, list):
        rem = grouped_points.copy()
        perms.append(find_subpath(None, rem))
    else:
        rem = grouped_points.copy()
        rem.remove(item)
        perms.append(find_path(item, rem))

# find shortest permutation
solution = []
minSolution = []
minDist = float('inf')
for perm in perms:
    dist = util.calcDist(perm)
    if dist < minDist:
        minSolution = perm
        minDist = dist

timeEnd = time.perf_counter()

for i in range(len(minSolution) - 1):
    wsp.ax[1].plot([minSolution[i].x, minSolution[i+1].x],[minSolution[i].y, minSolution[i+1].y], color="red")
wsp.ax[0].set_title(f"#WSP={wsp_count}, s={s}")
wsp.ax[1].set_title(f"TSP Path: n={len(points)}, length={minDist:0.4f}")

print("")
print("Solution:", minSolution)
print("Solution Distance:", minDist)
print(len(perms), "permutations examined")
print(f"Solution found in {timeEnd - timeStart:0.4f} seconds")'''
print("___________________________________________________________")
plt.show()