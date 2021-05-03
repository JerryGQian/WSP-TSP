from wsp import wsp_prq
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import time

# run algorithm
# >> python tsp-nnp.py <points file> <separation factor> <flags:{-d, -bf}>
# >> python tsp-nnp.py att48.tsp 1 -d
# -d: debug info quadtree and WSP
# -bf: brute force (turns off WSP)

timeInit = time.perf_counter()

filename = "data/points.txt"
s = 1           # default separation factor
wsp_mode = True # uses WSPs
debug = False   # debug info for Quadtree and WSPs

if (len(sys.argv) >= 2):
  filename = "data/" + sys.argv[1]
if (len(sys.argv) >= 3):
  s = int(sys.argv[2])
# check flags
for arg in sys.argv:
    if arg == "-d":
        debug = True
    if arg == "-bf":
        wsp_mode = False
# build WSP tree
wspTreeNode = wsp_prq.runWSP(filename, s, debug)

timeStart = time.perf_counter()

# calculate well separated dictionary
ws = dict() # point -> set of well separated points (far away by WSP)
ws_orig = dict() # point a -> dict( WSP point b -> WSP set containing point a )
points = wspTreeNode.get_points()
wsp_count = 0
for p in points:
    ws[p] = set()
    ws_orig[p] = dict()

q = [wspTreeNode]
if wsp_mode:
    while len(q) > 0:
        anode = q[0]
        q = q[1:]
        # add WSP relationships to ws
        if anode.connection != None:
            wsp_count += 1
            bnode = anode.connection
            apoints = anode.get_points()
            bpoints = bnode.get_points()
            for a in apoints:
                for b in bpoints:
                    ws[a].add(b)
                    ws[b].add(a)
                    ws_orig[a][b] = apoints
                    ws_orig[b][a] = bpoints
        if anode.divided:
            q.append(anode.ne)
            q.append(anode.nw)
            q.append(anode.sw)
            q.append(anode.se)
            
points = wspTreeNode.get_points()
num_points = len(points)
print("___________________________________________________________")
print(num_points, "points")
print(int(wsp_count/2), "WSPs found")
print(f"Loaded in {timeStart - timeInit:0.4f} seconds")
#print(ws)

# traversal
solution = []
minSolution = []
minDist = float('inf')


def euclidDist(p1, p2):
    return math.sqrt( ((p2.x - p1.x) ** 2) + ((p2.y - p1.y) ** 2) )

def calcDist(points):
    dist = 0
    for i in range(len(points) - 1):
        dist += euclidDist(points[i], points[i+1])
    return dist

def findPath(start, rem):
    perm = [start]
    if len(perm) == num_points:
        return [perm]
    while len(rem) > 0:
        minNext = None
        minNextDist = float('inf')
        for r in rem:
            last_point = perm[len(perm) - 1]
            orig_set_finished = True
            if r in ws[last_point]: # checks if all points in last_point <-> r set have been visited
                for p in perm:
                    if (r in ws_orig[last_point]) and (p in ws_orig[last_point][r]):
                        orig_set_finished = False
            if (r not in ws[last_point]) or (orig_set_finished) or not wsp_mode:
                curNextDist = euclidDist(last_point, r)
                if curNextDist < minNextDist:
                    minNext = r
                    minNextDist = curNextDist
        if minNext == None:
            minNext = rem[0]
        #print(minNext, rem)
        perm.append(minNext)
        rem.remove(minNext)
    return perm

perms = []
for p in points:
    rem = points.copy()
    rem.remove(p)
    perms.append(findPath(p, rem))


for perm in perms:
    dist = calcDist(perm)
    if dist < minDist:
        minSolution = perm
        minDist = dist

timeEnd = time.perf_counter()
print("")
print("Solution:", minSolution)
print("Solution Distance:", minDist)
print(len(perms), "permutations examined")
print(f"Solution found in {timeEnd - timeStart:0.4f} seconds")
print("___________________________________________________________")
plt.show()
#print(calcDist([wsp.Point(1,0),wsp.Point(8,1)]))