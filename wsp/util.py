import math

def euclidDist(p1, p2):
    return math.sqrt(((p2.x - p1.x) ** 2) + ((p2.y - p1.y) ** 2))

def calcDist(points):
    dist = 0
    for i in range(len(points) - 1):
        dist += euclidDist(points[i], points[i+1])
    return dist