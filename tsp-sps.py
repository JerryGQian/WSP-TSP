from wsp import wsp
from wsp import ds
from wsp import util
from wsp import cmd_parse
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import time

# run algorithm
# >> python tsp-sps.py <points file> <separation factor> <quadtree:{-pr, -point/-p}> <flags:{-d, -bf}>
# >> python tsp-sps.py att48.tsp 1 -pr
# -d: debug info quadtree and WSP
# -bf: brute force (turns off WSP)

timeInit = time.perf_counter()

filename, s, wsp_mode, debug, quadtree, bucket = cmd_parse.parse_cmd(sys.argv)
# build WSP tree
wspTreeNode, wsp_count = wsp.runWSP(filename, s, debug, quadtree, bucket)

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

def find_relations(tree_node, add=True):
    sub_relations = set()

    if len(tree_node.connection) > 0:
        for node in tree_node.connection:
            for p in node.get_points():
                sub_relations.add(p)

    if tree_node.divided:
        sub_relations.union(find_relations(tree_node.ne, True))
        sub_relations.union(find_relations(tree_node.nw, True))
        sub_relations.union(find_relations(tree_node.sw, True))
        sub_relations.union(find_relations(tree_node.se, True))

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
        if add:
            #print("insert:", (node_point_set, sub_relations.copy()))
            splits.insert(0, (node_point_set, sub_relations.copy()))

    '''if tree_node.divided:
        find_relations(tree_node.ne, True)
        find_relations(tree_node.nw, True)
        find_relations(tree_node.sw, True)
        find_relations(tree_node.se, True)'''

    return sub_relations

find_relations(wspTreeNode, True)
print(splits)

def apply_split(pair, glist):
    list1 = []
    list2 = []    
    to_remove = []
    to_add = []
    for item in glist:
        if isinstance(item, list):
            # recurse down list, add later
            to_remove.append(item)
            to_add.append(apply_split(pair, item.copy()))
        else:
            if item in pair[0]:
                list1.append(item)
                #glist.remove(item)
            elif item in pair[1]:
                list2.append(item)
                #glist.remove(item)
    for item in (list1 + list2 + to_remove):
        glist.remove(item)
    for item in to_add: # sublists
        glist.append(item)

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
    if len(pair[0]) >= 2 or len(pair[1]) >= 2:
        grouped_points = apply_split(pair, grouped_points.copy())
        #print("after pair:", pair, " -> ", grouped_points)
    #print("GP: ", grouped_points)
#print(points)
#print(grouped_points)

# traversal
'''def closest_sub(p1, sublist):
  """Min dist point of sublist from p1"""
  mind = 99999999
  minPoint = None
  for p2 in sublist:
    dist = p1.distance_to(p2)
    if dist < mind:
        mind = dist
        minPoint = p2
  return minPoint, mind


def find_pathX(start, rem):
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
                curNextDist = last_point.distance_to(r)
            if curNextDist < minNextDist:
                minNext = r
                minNextDist = curNextDist
        if minNext == None:
            minNext = rem[0] # shouldnt be needed?
        if isinstance(minNext, list):
            #print("minNext sub", minNext)
            path += find_subpath(last_point, minNext)
            #print(minNext)
        #else:
            #print("should be point",minNext, len(rem))
            
            path.append(minNext)
        rem.remove(minNext)
    return path'''

def does_list_contain(point, lst):
    for item in lst:
        if isinstance(item, list):
            if does_list_contain(point, item):
                return True
        else:
            if item == point:
                return True
    return False

def find_list_with(point, lst):
    for item in lst:
        if isinstance(item, list):
            if does_list_contain(point, item):
                return item
        else:
            if item == point:
                return item
    return None

def clean_deep_lists(lst):
    if len(lst) == 1:
        lst[0] = clean_deep_lists(lst[0])
        return lst[0]
    return lst

def find_path(start, glist, end):

    #print(start, end)
    reverse = False
    if start == None:
        reverse = True
        start = end
        end = None
        #rem = glist
    orig_start = start
    glist = clean_deep_lists(glist)

    #print("glist:", glist)
    #print("start:", start)
    rem = glist.copy()
    if start not in glist:
        start = find_list_with(start, glist)
        #print("new start", start)
    rem.remove(start) 
    #print("rem:  ", rem)
    #print(len(rem))
    trio_path = [(None, start, None)]
    if len(rem) == 0:
        trio_path = [(orig_start, start, None)]
    # perform nearest neighbor on topmost layer
    while len(rem) > 0:
        prev_trio = trio_path[len(trio_path) - 1] # trio is (prev point, item, next point)
        prev_item = prev_trio[1]
        minNext = None
        minNextDist = float('inf')
        for r in rem: # look through rem for next
            p_A, p_B = util.min_proj_set_or_point(prev_item, r)
            dist = p_A.distance_to(p_B)
            if end != None and r == end:
                dist += 9999999999 # end_penalty
            if dist < minNextDist:
                minNext = (p_A, r, p_B) # (prev point, next item, next point)
                minNextDist = dist
        #print("conn next", minNext)
        #print("triopath", trio_path, minNext, rem, end)
        trio_path[len(trio_path) - 1] = (prev_trio[0], prev_trio[1], minNext[0])
        trio_path.append((minNext[2], minNext[1], None))
        rem.remove(minNext[1])

    # convert to list, make subcalls
    path = []
    print("trio_path:", trio_path)
    for trio in trio_path:
        #print("trio:", trio)
        if isinstance(trio[1], list):
            path += find_path(trio[0], trio[1], trio[2])
        else:
            path.append(trio[1])
    if reverse:
        path.reverse()
    return path


# search for permutations
perms = []
for item in grouped_points:
    rem = grouped_points.copy()
    #rem.remove(item)
    print("start", item)
    perms.append(find_path(item, rem, item))
    '''if isinstance(item, list):
    else:
        perms.append(find_path(item, rem))'''
'''for item in grouped_points:
    if isinstance(item, list):
        rem = grouped_points.copy()
        perms.append(find_subpath(None, rem))
    else:
        rem = grouped_points.copy()
        rem.remove(item)
        perms.append(find_path(item, rem))'''

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