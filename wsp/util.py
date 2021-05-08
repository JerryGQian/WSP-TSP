import math

def euclidDist(p1, p2):
    return math.sqrt(((p2.x - p1.x) ** 2) + ((p2.y - p1.y) ** 2))

def calcDist(points):
    dist = 0
    for i in range(len(points) - 1):
        dist += euclidDist(points[i], points[i+1])
    return dist

def sublist_get_points(lst):
    points = []
    for item in lst:
        if isinstance(item, list):
            points += sublist_get_points(item)
        else:
            points.append(item)
    return points

def min_proj_set_or_point(item_A, item_B):
    if not isinstance(item_A, list):
        item_A = [item_A]
    if not isinstance(item_B, list):
        item_B = [item_B]
    item_A = sublist_get_points(item_A)
    item_B = sublist_get_points(item_B)
    return min_proj(item_A, item_B)

def min_proj(set_A, set_B):
  """Min pair between points from set_A and set_B"""
  mind = 99999999
  min_p1 = None
  min_p2 = None
  for p_A in set_A:
    for p_B in set_B:
      dist = p_A.distance_to(p_B)
      if dist < mind:
        mind = dist
        min_p1 = p_A
        min_p2 = p_B
  return min_p1, min_p2