from wsp import ds
import sys


def parse_cmd(argv):
    filename = "data/custom1.txt"
    s = 1           # default separation factor
    wsp_mode = True # uses WSPs
    debug = False   # debug info for Quadtree and WSPs
    quadtree = ds.PMRQuadTree
    bucket = 1

    if (len(argv) >= 2):
        filename = "data/" + argv[1]
    if (len(argv) >= 3):
        s = float(argv[2])
    # check flags
    for arg in argv:
        if arg == "-pkpr" or arg == "-pk":
            quadtree = ds.PKPRQuadTree
        if arg == "-pmr":
            quadtree = ds.PMRQuadTree
        if arg == "-pr":
            quadtree = ds.PRQuadTree
        if arg == "-point" or arg == "-p":
            quadtree = ds.PointQuadTree
        if arg[:2] == "-b":
            bucket = int(arg[2:])
        if arg == "-d":
            debug = True
    
    return filename, s, wsp_mode, debug, quadtree, bucket

