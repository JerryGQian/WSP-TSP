from wsp import wsp
from wsp import ds
from wsp import util
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import time

# run algorithm
# >> python tsp-spp.py <points file> <separation factor> <quadtree:{-pr, -point/-p}> <flags:{-d, -bf}>
# >> python tsp-spp.py att48.tsp 1 -p -d
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
wspTreeNode, wsp_count = wsp.runWSP(filename, s, debug, quadtree)

timeStart = time.perf_counter()

# calculate well separated dictionary
ws = dict() # point -> set of well separated points (far away by WSP)
ws_orig = dict() # point a -> dict( WSP point b -> WSP set containing point a )
points = wspTreeNode.get_points()
for p in points:
    ws[p] = set()
    ws_orig[p] = dict()

points = wspTreeNode.get_points()
num_points = len(points)
print("___________________________________________________________")
print(num_points, "points")
print(int(wsp_count), "WSPs found")
print(f"Loaded in {timeStart - timeInit:0.4f} seconds")

# [ ([a,b,c],[d,e,f]) , ([a,b,c],[d,e,f]) ]
splits = []

def find_relations(tree_node):
    sub_relations = set()

    if len(tree_node.connection) > 0:
        for node in tree_node.connection:
            for p in node.get_points():
                sub_relations.add(p)

    if tree_node.divided:
        sub_relations.union(find_relations(tree_node.ne))
        sub_relations.union(find_relations(tree_node.nw))
        sub_relations.union(find_relations(tree_node.sw))
        sub_relations.union(find_relations(tree_node.se))

    node_point_set = set(tree_node.get_points())
    to_remove = []
    if len(sub_relations) > 0:
        for p in sub_relations:
            if p in node_point_set:
                #print("removing")
                to_remove.append(p)
                #sub_relations.remove(p)
        for p in to_remove:
            sub_relations.remove(p)
        splits.insert(0, (node_point_set, sub_relations.copy()))

    return sub_relations

find_relations(wspTreeNode)
print(splits)

def apply_split(pair, glist):
    list1 = []
    list2 = []    
    for item in glist:
        if isinstance(item, list):
            item = apply_split(pair, item.copy())
        else:
            if item in pair[0]:
                list1.append(item)
                #glist.remove(item)
            elif item in pair[1]:
                list2.append(item)
                #glist.remove(item)
    for item in (list1 + list2):
        glist.remove(item)

    if len(list1) == 1:
        glist.append(list1[0])
    elif len(list1) > 0:
        glist.append(list1)
    if len(list2) == 1:
        glist.append(list2[0])
    elif len(list2) > 0:
        glist.append(list2)
    return glist

grouped_points = points.copy()

for pair in splits:
    grouped_points = apply_split(pair, grouped_points.copy())
    #print("GP: ", grouped_points)
#print(points)
print(grouped_points)
'''for g in grouped_points:
    if isinstance(g, list):
        for i in g:
            if isinstance(i, list):
                print("has sublist")
        print(len(g))
    else:
        print("1")'''

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
            print("minNext sub", minNext)
            path += find_subpath(last_point, minNext)
            #print(minNext)
        else:
            print("should be point",minNext, len(rem))
            
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
print("perms", perms)
timeEnd = time.perf_counter()

for i in range(len(minSolution) - 1):
    wsp.ax[1].plot([minSolution[i].x, minSolution[i+1].x],[minSolution[i].y, minSolution[i+1].y], color="red")
wsp.ax[0].set_title(f"#WSP={wsp_count}, s={s}")
wsp.ax[1].set_title(f"TSP Path: n={len(points)}, length={minDist:0.4f}")

print("")
print("Solution:", minSolution)
print("Solution Distance:", minDist)
print(len(perms), "permutations examined")
print(f"Solution found in {timeEnd - timeStart:0.4f} seconds")
print("___________________________________________________________")
plt.show()