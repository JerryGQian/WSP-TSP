from wsp import file_load
from wsp import ds
import numpy as np
import sys
import matplotlib.pyplot as plt

# USES PR QUADTREE!
fig, ax = plt.subplots(1, 2, figsize=(12,6))

def runWSP(filename, s, debug, quadtree, bucket):
  points, minX, minY, maxX, maxY = file_load.loadFromFile(filename)
  # build point quadtree, insert in order
  rootNode = quadtree(ds.Rect(minX,minY,maxX,maxY), ax, bucket)

  #rootNode = QuadTree(Rect(bounds[0],bounds[1],bounds[2],bounds[3]))
  for point in points:
    rootNode.insert(point)

  if (debug):
    print(points, "\n")
    print(rootNode, "\n\n")

  wsp_count = 0

  # WSP queue search
  #s = 1 # wsp separation factor
  ax[0].plot([rootNode.boundary.xMin, rootNode.boundary.xMax],[rootNode.boundary.yMin, rootNode.boundary.yMin], color="gray")
  ax[0].plot([rootNode.boundary.xMin, rootNode.boundary.xMax],[rootNode.boundary.yMax, rootNode.boundary.yMax], color="gray")
  ax[0].plot([rootNode.boundary.xMin, rootNode.boundary.xMin],[rootNode.boundary.yMin, rootNode.boundary.yMax], color="gray")
  ax[0].plot([rootNode.boundary.xMax, rootNode.boundary.xMax],[rootNode.boundary.yMin, rootNode.boundary.yMax], color="gray")
  ax[1].plot([rootNode.boundary.xMin, rootNode.boundary.xMax],[rootNode.boundary.yMin, rootNode.boundary.yMin], color="gray")
  ax[1].plot([rootNode.boundary.xMin, rootNode.boundary.xMax],[rootNode.boundary.yMax, rootNode.boundary.yMax], color="gray")
  ax[1].plot([rootNode.boundary.xMin, rootNode.boundary.xMin],[rootNode.boundary.yMin, rootNode.boundary.yMax], color="gray")
  ax[1].plot([rootNode.boundary.xMax, rootNode.boundary.xMax],[rootNode.boundary.yMin, rootNode.boundary.yMax], color="gray")
  queue = [(rootNode, rootNode)]
  while len(queue) > 0:
    pair = queue[0]
    queue = queue[1:]

    block_A, block_B = pair[0], pair[1]

    if len(block_A) == 0 or len(block_B) == 0:
      continue
    
    points_A = block_A.get_points()
    points_B = block_B.get_points()
    if ds.min_dist(block_A, block_B) >= s * block_A.diameter() or (len(points_A) == 1 and len(points_B) == 1 and not block_A.divided  and not block_B.divided):
      if (debug):
        print("found a WSP: ", block_A.str_short(), " <~~~~~> ", block_B.str_short())
      wsp_count += 1
      block_A.connection.append(block_B)
      block_B.connection.append(block_A)
      circle1 = plt.Circle(block_A.center(), block_A.diameter() / 2, color='r', fill=False)
      circle2 = plt.Circle(block_B.center(), block_B.diameter() / 2, color='r', fill=False)
      ax[0].add_patch(circle1)
      ax[0].add_patch(circle2)
      #line
      ax[0].plot([block_A.center()[0], block_B.center()[0]],[block_A.center()[1], block_B.center()[1]])
      continue

    if block_A.diameter() > block_B.diameter():
      if (block_A.divided):
        queue.append((block_B, block_A.nw))
        queue.append((block_B, block_A.ne))
        queue.append((block_B, block_A.se))
        queue.append((block_B, block_A.sw))
    else:
      if (block_B.divided):
        queue.append((block_A, block_B.nw))
        queue.append((block_A, block_B.ne))
        queue.append((block_A, block_B.se))
        queue.append((block_A, block_B.sw))
  
  # plot points
  x = []
  y = []
  for p in points:
    x.append(p.x)
    y.append(p.y)
  #fig = plt.scatter(x, y)
  ax[0].scatter(x, y)
  ax[1].scatter(x, y)

  return rootNode, wsp_count