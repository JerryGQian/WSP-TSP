from wsp import ds
import numpy as np
import random

def loadFromFile(filename):
    points = []
    # read points from file
    bounds = []
    with open(filename, 'r') as f: # reads .TSP files
        line = f.readline()
        mode = "start"
        while line != '':  # The EOF char is an empty string
            if line == "EOF\n":
                break
            if mode == "start":
                if line[len(line) - 1] == "\n":
                    line = line[:-1]
                if len(line) > 7 and line[:7] == "bounds:":
                    bounds = [int(i) for i in line[7:].split()]
                if line == "NODE_COORD_SECTION":
                    mode = "node"
                #if len(line) == 0 or line[0] == '#': # ignores empty lines and #comments
                line = f.readline()
                continue
            # start reading node coords
            if mode == "node":
                if line == "TOUR_SECTION\n":
                    mode = "tour"
                    break
                splitLine = line.split()
                if len(splitLine) == 3:
                    splitLine = splitLine[1:]
                p = ds.Point(float(splitLine[0].strip()), float(splitLine[1].strip()))
                points.append(p)
                line = f.readline()

    # find boundaries
    minX = float('inf')
    minY = float('inf')
    maxX = float('-inf')
    maxY = float('-inf')
    for p in points:
        if p.x < minX:
            minX = p.x
        if p.y < minY:
            minY = p.y
        if p.x > maxX:
            maxX = p.x
        if p.y > maxY:
            maxY = p.y
    maxX += 1
    maxY += 1

    random.shuffle(points)

    return (points, minX, minY, maxX, maxY)