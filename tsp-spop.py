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

filename, s, wsp_mode, debug, shrink, quadtree, bucket = cmd_parse.parse_cmd(sys.argv)
# build WSP tree
wspTreeNode, wsp_count = wsp.runWSP(filename, s, debug, shrink, quadtree, bucket)

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

    if quadtree == ds.PKPRQuadTree or quadtree == ds.PKPMRQuadTree:
        for child in tree_node.children:
            sub_relations.union(find_relations(child, True))
    else:
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
                to_remove.append(p)
        for p in to_remove:
            sub_relations.remove(p)
        if add:
            if len(splits) == 0:
                splits.append((node_point_set, sub_relations.copy()))
            else:
                for i in range(len(splits)):
                    #print(len(node_point_set), len(splits[i][0]))
                    len_score1 = (len(node_point_set) + len(sub_relations)) / 2 - abs(len(node_point_set) - len(sub_relations)) / .5
                    len_score2 = (len(splits[i][0]) + len(splits[i][0])) / 2 - abs(len(splits[i][0]) - len(splits[i][0])) / .5
                    if len_score1 > len_score2:
                        splits.insert(i, (node_point_set, sub_relations.copy()))
                        break

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
#print(points)
print("grouped_points", grouped_points)

# traversal

def does_list_contain(point, lst):
    for item in lst:
        if isinstance(item, list):
            if does_list_contain(point, item):
                return True
        else:
            if item == point:
                return True
    return False

def get_points_from_list(lst):
    return get_points_from_list_rec(lst, [])

def get_points_from_list_rec(lst, points):
    for item in lst:
        if isinstance(item, list):
            get_points_from_list_rec(item, points)
        else:
            points.append(item)
    return points

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
    if isinstance(lst, list) and len(lst) == 1:
        lst[0] = clean_deep_lists(lst[0])
        return lst[0]
    return lst

def find_path(start, glist, end, depth=-1, init=False):
    #if depth > 10:
    #    return None
    reverse = False
    if start == None:
        reverse = True
        start = end
        end = None
    #orig_start = start
    glist = clean_deep_lists(glist)

    start_item = start
    end_item = end
    if init:
        #start_item = start
        #end_item = end
        start = None
        end = None
        end_item = None
        #start, end = util.min_proj_set_or_point(start, end)
    else:
        if start not in glist:
            start_item = find_list_with(start, glist)
        if end not in glist:
            end_item = find_list_with(end, glist)

    print("start", start, "end", end)
    #print("start_item", start_item, "end_item", end_item)
    #print("glist", glist)
    
    rem = glist.copy()
    #print(isinstance(start_item, list) and end != None and does_list_contain(end, start_item))
    if not init and ((isinstance(start_item, list) and end != None and does_list_contain(end, start_item)) or (end_item == start_item)):
        # if start and end in same sublist, disassemble start sublist into rem
        print("BREAKING up dual connection subproblem")
        rem.remove(start_item)
        subitems = get_points_from_list(start_item)

        if start in subitems:
            start_to_add = start
        else:
            start_to_add = find_list_with(start, subitems)

        if end in subitems:
            end_to_add = end
        else:
            end_to_add = find_list_with(end, subitems)
        #print("subitems bef:", subitems)
        subitems.remove(start_to_add)
        subitems.remove(end_to_add)

        #print("starttoadd", start_to_add, "endtoadd", end_to_add, "subitems", subitems)
        #print("subitems after:", subitems)
        if len(subitems) > 0:
            rem.append(subitems)#get_points_from_list(start_item)
        #rem.remove(start)
        #rem.append(start_to_add)
        rem.append(end_to_add)

        start_item = start_to_add
        end_item = end_to_add
    else:    
        # regular remove behavior
        rem.remove(start_item) 
    print("rem", rem)
    trio_path = [(start, start_item, None)]

    #if len(rem) == 0:
    #    #print("at the end, trio:",  start)
    #    trio_path = [(orig_start, start, None)]
    def point_list_contains(pointlist, point):
        if isinstance(pointlist, list):
            return does_list_contain(point, pointlist)
        else:
            return point == pointlist

    # perform nearest neighbor on topmost layer
    #min_proj_dict = {}
    def buildPerms(perm, rem):
        #print(perm, rem)
        next_perms = []
        if len(rem) == 0:
            return [perm]
        #print(len(rem))
        for r in rem:
            last_point = perm[len(perm) - 1]
            new_point_list = perm.copy()
            new_point_list.append(r)
            new_rem = rem.copy()
            new_rem.remove(r)
            if len(rem) == 0:
                next_perms.append(new_point_list)
            else:
                next_perms.extend(buildPerms(new_point_list, new_rem))
        return next_perms

    # start at points[0]
    #print("Building perms!", len(rem))
    if len(rem) < 9:
        print("BRUTE FORCE depth:", depth)
        perms = buildPerms([start_item], rem)
        #print("Perms:", len(perms))
        best_perm = None
        min_d = float('inf')
        for perm in perms:
            #if end != None and end == orig_start:
                #perm.append(end)
                #print("appended perm", perm)
            #print("Perm:",perm)
            #print("compare: ", start, perm[0], end, perm[len(perm) - 1], len(perms))
            #print("inlistL ", does_list_contain(orig_start, perm[0]), does_list_contain(end, perm[len(perm) - 1]))
            #print(does_list_contain(orig_start, perm[0]), does_list_contain(end, perm[len(perm) - 1]))
            test_perm = init
            if not test_perm:
                start_cond = False
                end_cond = False
                if isinstance(perm[0], list):
                    start_cond = does_list_contain(start, perm[0])
                else:
                    start_cond = start == perm[0]

                if isinstance(perm[len(perm) - 1], list):
                    end_cond = does_list_contain(end, perm[len(perm) - 1]) or end == None
                else:
                    end_cond = end == perm[len(perm) - 1] or end == None

                test_perm = start_cond and end_cond
            if init or test_perm:
                #print("not skip")
                dist = 0
                for i in range(len(perm) - 1):
                    #print("TESTING", perm[i], perm[i+1])
                    p_A, p_B = util.min_proj_set_or_point(perm[i], perm[i+1])
                    dist += p_A.distance_to(p_B)
                #print(perm, dist)
                if dist < min_d:
                    min_d = dist
                    best_perm = perm
            #else:
                #print("skipping")
        
        #best_perm
        #print("best perm - depth:", depth, min_d, best_perm)

        #trio_path = []
        for i in range(1,len(best_perm)):
            item = best_perm[i]
            #print("item in bestperm", item)
            #if len(trio_path) > 0:
            prev_trio = trio_path[len(trio_path) - 1]
            prev_item = prev_trio[1]
            p_prev, p_cur, p2_prev, p2_cur = util.min_proj_set_or_point(prev_item, item, 2)
            if p_cur == end and i == len(best_perm) - 1:
                p_cur = p2_cur
            trio_path[len(trio_path) - 1] = (prev_trio[0], prev_trio[1], p_prev)
            item = clean_deep_lists(item)
            trio_path.append((p_cur, item, None))
            #print((p_cur, item, None))
            #print("trio_path 1:", trio_path)
            #else:
            #    item = clean_deep_lists(item)
            #    trio_path.append((None, item, None)) ####

        #if best_perm[0] == best_perm[len(best_perm) - 1]:
            #trio_path[0] = ()
        #print("bf triopath", trio_path)
    else:
        print("NEAREST NEIGHBOR depth:", depth, len(rem))
        while len(rem) > 0:
            prev_trio = trio_path[len(trio_path) - 1] # trio is (prev point, item, next point)
            prev_item = prev_trio[1]
            minNext = None
            minNextDist = float('inf')
            #print("END", end, "END ITEM:",  end_item, "prev", prev_item)
            #print("rem", rem)
            for r in rem: # look through rem for next
                p_A, p_B = util.min_proj_set_or_point(prev_item, r)

                dist = p_A.distance_to(p_B)
                if dist < minNextDist and (len(rem) == 1 or r != end_item):
                    #print("checking", r, len(rem), point_list_contains(r, end_item))
                    minNext = (p_A, r, p_B) # (prev point, next item, next point)
                    minNextDist = dist
                #else:
                    #print("skipping", r, len(rem), point_list_contains(r, end_item))
            #print("minNext", minNext)
            trio_path[len(trio_path) - 1] = (prev_trio[0], prev_trio[1], minNext[0])
            trio_path.append((minNext[2], minNext[1], None))
            rem.remove(minNext[1])

    # assign end to end
    if init:
        first_trio = trio_path[0]
        last_trio = trio_path[len(trio_path) - 1]
        p_start, p_end, p_start2, p_end2 = util.min_proj_set_or_point(start_item, last_trio[1], 2)
        if first_trio[2] == p_start:
            p_start = p_start2
        if last_trio[0] == p_end:
            p_end = p_end2
        trio_path[0] = (p_start, first_trio[1], first_trio[2])
        trio_path[len(trio_path) - 1] = (last_trio[0], last_trio[1], p_end)
    elif end != None: # ?init
        #print("assigning end", end)
        last_trio = trio_path[len(trio_path) - 1]
        trio_path[len(trio_path) - 1] = (last_trio[0], last_trio[1], end)

    # reroute trios with same start and end to second min projection
    for i in range(len(trio_path)-1):
        trio = trio_path[i]
        if trio[0] == trio[2]: # fix if same
            next_trio = trio_path[i + 1]
            #print("fixing...", trio[1], next_trio[1])
            _, _, p_A, _ = util.min_proj_set_or_point(trio[1], next_trio[1], True)
            #print("pA", p_A)
            if p_A != None:
                print("fixing trios", trio[0], trio[2], "=>", p_A)
                trio_path[i] = (trio[0], trio[1], p_A)
    for i in range(len(trio_path)):
        trio = trio_path[i]
        if not isinstance(trio[1], list):
            trio_path[i] = (trio[1], trio[1], trio[1])

    print("trio path:", trio_path)
    # convert to list, make subcalls
    path = []
    last_trio = trio_path[len(trio_path) - 1]
    #trio_path[len(trio_path) - 1] = (last_trio[0], last_trio[1], end) ## why this here?
    #if depth < 2 or True:
    #    print(depth, "trio_path:", trio_path)
    for trio in trio_path:
        #print("trio:", trio)
        if isinstance(trio[1], list):
            #print("testing trios", trio[0], trio[1], trio[2])
            npath = find_path(trio[0], trio[1], trio[2], depth + 1)
            #if depth < 2 or True:
            #    #print("calling", trio)
            #    print(depth + 1, "start and ends:", trio[0], trio[2])
            #    print(depth + 1, "res path:", reverse, npath)
            print(npath, trio[0], trio[1], trio[2])
            path += npath
        else:
            #print(depth + 1, "res point:", trio[1])
            path.append(trio[1])
    if reverse:
        path.reverse()
    return path


# search for permutations
rem = grouped_points.copy()
print(grouped_points[0])
perm = find_path(grouped_points[0], rem, grouped_points[0], 0, True)
perm.append(perm[0])
print(perm)
'''perms = []
for item in grouped_points:
    rem = grouped_points.copy()
    perm = find_path(item, rem, item, 0)
    perm.append(perm[0])
    perms.append(perm)'''

# find shortest permutation
solution = []
minSolution = []
minDist = float('inf')
'''for perm in perms:
    print(perm)
    dist = util.calcDist(perm)
    if dist < minDist:
        minSolution = perm
        minDist = dist'''
#print("perms", perms)
minSolution = perm
minDist = util.calcDist(perm)
timeEnd = time.perf_counter()

for i in range(len(minSolution) - 1):
    wsp.ax[1].plot([minSolution[i].x, minSolution[i+1].x],[minSolution[i].y, minSolution[i+1].y], color="red")
wsp.ax[0].set_title(f"#WSP={wsp_count}, s={s}")
wsp.ax[1].set_title(f"TSP Path: n={len(points)}, length={minDist:0.4f}")

print("")
print("Solution:", minSolution)
print("Solution Distance:", minDist)
#print(len(perms), "permutations examined")
print(f"Solution found in {timeEnd - timeStart:0.4f} seconds")
print("___________________________________________________________")
plt.show()