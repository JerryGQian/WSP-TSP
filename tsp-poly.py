from wsp import wsp
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# run algorithm
# >> python tsp-poly.py <points file> <separation factor> <flags:{-d, -bf}>
# >> python tsp-poly.py points1.txt 1 -d
# -d: debug info quadtree and WSP
# -bf: brute force (turns off WSP)

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
wspTreeNode = wsp.runWSP(filename, s, debug)

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

def buildPerms(perm, rem):
    next_perms = []
    if len(perm) == num_points:
        return [perm]
    for r in rem:
        last_point = perm[len(perm) - 1]
        orig_set_finished = True
        if r in ws[last_point]: # checks if all points in last_point <-> r set have been visited
            for p in perm:
                if (r in ws_orig[last_point]) and (p in ws_orig[last_point][r]):
                    orig_set_finished = False
        if (r not in ws[last_point]) or (orig_set_finished):
            new_point_list = perm.copy()
            new_point_list.append(r)
            new_rem = rem.copy()
            new_rem.remove(r)
            if len(new_point_list) == num_points:
                next_perms.append(new_point_list)
            else:
                next_perms.extend(buildPerms(new_point_list, new_rem))
    return next_perms

perms = []
for p in points:
    rem = points.copy()
    rem.remove(p)
    perms.extend(buildPerms([p], rem))

print(len(perms), "permutations examined")

for perm in perms:
    dist = calcDist(perm)
    if dist < minDist:
        minSolution = perm
        minDist = dist
print("")
print("Solution:", minSolution)
print("Solution Distance:", minDist)
print("___________________________________________________________")
plt.show()
#print(calcDist([wsp.Point(1,0),wsp.Point(8,1)]))