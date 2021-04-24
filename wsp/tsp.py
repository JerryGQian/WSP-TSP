import wsp
import sys

# run algorithm
# >> python tsp.py points.txt 1
filename = "points.txt"
s = 1
debug = False
if (len(sys.argv) >= 2):
  filename = sys.argv[1]
if (len(sys.argv) >= 3):
  s = int(sys.argv[2])
if (len(sys.argv) >= 4):
    if sys.argv[3] == "-d":
        debug = True
# build WSP tree
wspTreeNode = wsp.runWSP(filename, s, debug)

# calculate well separated dictionary
ws = dict()
points = wspTreeNode.get_points()
for p in points:
    ws[p] = set()
q = [wspTreeNode]
while len(q) > 0:
    anode = q[0]
    q = q[1:]
    # add WSP relationships to ws
    if anode.connection != None:
        bnode = anode.connection
        apoints = anode.get_points()
        bpoints = bnode.get_points()
        for a in apoints:
            for b in bpoints:
                ws[a].add(b)
                ws[b].add(a)
    if anode.divided:
        q.append(anode.ne)
        q.append(anode.nw)
        q.append(anode.sw)
        q.append(anode.se)

print(ws)

# traversal
solution = []
minSolution = []
minDist = float('inf')

points = wspTreeNode.get_points()
for p in points:
    x = 0

print("Solution:", minSolution)
print("Solution Distance:", minDist)