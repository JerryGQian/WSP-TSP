import numpy as np

class Point:
    #A point located at (x,y) in 2D space.
    def __init__(self, x, y):
        self.x, self.y = x, y
    def __repr__(self):
        return '{}'.format(str((self.x, self.y)))
    def __str__(self):
        return 'P({:.2f}, {:.2f})'.format(self.x, self.y)
    def __add__(self, o):
        return Point(self.x + o.x, self.y + o.y)
    def __sub__(self, o):
        return Point(self.x - o.x, self.y - o.y)
    def __truediv__(self, o):
        return Point(self.x / o, self.y / o)
    def to_list(self):
        return [self.x, self.y]
    def distance_to(self, other):
        try:
            other_x, other_y = other.x, other.y
        except AttributeError:
            other_x, other_y = other
        return np.hypot(self.x - other_x, self.y - other_y)

class Rect:
    def __init__(self, xMin, yMin, xMax, yMax):
      self.xMin = xMin
      self.yMin = yMin
      self.xMax = xMax
      self.yMax = yMax

    def __repr__(self):
        return str((self.xMin,self.yMin, self.xMax, self.yMax))

    def __str__(self):
        return '({:.2f}, {:.2f}, {:.2f}, {:.2f})'.format(self.xMin,
                    self.yMin, self.xMax, self.yMax)

    def contains(self, point):
        """Is point (a Point object or (x,y) tuple) inside this Rect?"""
        try:
            point_x, point_y = point.x, point.y
        except AttributeError:
            point_x, point_y = point

        return (point_x >= self.xMin and
                point_x <  self.xMax and
                point_y >= self.yMin and
                point_y < self.yMax)

    def diameter(self):
      # diagonal
      return np.hypot(self.xMax - self.xMin, self.yMax - self.yMin)
    def center(self):
      return ((self.xMax + self.xMin) / 2, (self.yMax + self.yMin) / 2)

def min_dist(block_A, block_B):
    min_p1, min_p2 = min_proj(block_A, block_B)
    return min_p1.distance_to(min_p2)

def min_proj(block_A, block_B):
  """Min dist between points from Quadtree block_A and Quadtree block_B"""
  points_A = block_A.get_points()
  points_B = block_B.get_points()
  mind = 99999999
  min_p1 = None
  min_p2 = None
  for p_A in points_A:
    for p_B in points_B:
      dist = p_A.distance_to(p_B)
      if dist < mind:
        mind = dist
        min_p1 = p_A
        min_p2 = p_B
  return min_p1, min_p2

# POINT REGION QUADTREE

class PRQuadTree:
    """Point Region Quadtree implementation."""

    def __init__(self, boundary, ax, depth=0):
        """Initialize this node of the quadtree."""
        self.boundary = boundary # boundary of current block
        self.ax = ax
        self.depth = depth # mostly for string visualization spacing
        self.point = None # center point
        self.connection = [] # WSP connection
        self.divided = False # flag for if divided into 4 child quads

    def __str__(self):
        """Return a string representation of this node, suitably formatted."""
        sp = ' ' * self.depth * 2
        s = str(self.boundary) + ' --> ' + str(self.point) 
        if not self.divided:
            return s
        return s + '\n' + '\n'.join([
                sp + 'nw: ' + str(self.nw), sp + 'ne: ' + str(self.ne),
                sp + 'se: ' + str(self.se), sp + 'sw: ' + str(self.sw)])

    def str_short(self):
      return str(self.boundary)

    def diameter(self):
      return self.boundary.diameter()
    def center(self):
      return self.boundary.center()

    def divide(self):
        """Divide (branch) this node by spawning four children nodes around a point."""
        mid = Point((self.boundary.xMin + self.boundary.xMax) / 2, (self.boundary.yMin + self.boundary.yMax) / 2)
        self.nw = PRQuadTree(Rect(self.boundary.xMin, mid.y, mid.x, self.boundary.yMax), self.ax, self.depth + 1)
        self.ne = PRQuadTree(Rect(mid.x, mid.y, self.boundary.xMax, self.boundary.yMax), self.ax, self.depth + 1)
        self.se = PRQuadTree(Rect(mid.x, self.boundary.yMin, self.boundary.xMax, mid.y), self.ax, self.depth + 1)
        self.sw = PRQuadTree(Rect(self.boundary.xMin, self.boundary.yMin, mid.x, mid.y), self.ax, self.depth + 1)
        self.divided = True
        # reinsert point
        point_to_reinsert = self.point
        self.point = None
        if point_to_reinsert != None:
          self.insert(point_to_reinsert)
        # draw
        self.ax[0].plot([mid.x, mid.x],[self.boundary.yMin, self.boundary.yMax], color="gray")
        self.ax[0].plot([self.boundary.xMin, self.boundary.xMax],[mid.y, mid.y], color="gray")
        self.ax[1].plot([mid.x, mid.x],[self.boundary.yMin, self.boundary.yMax], color="gray")
        self.ax[1].plot([self.boundary.xMin, self.boundary.xMax],[mid.y, mid.y], color="gray")

    def insert(self, point):
        """Try to insert Point point into this QuadTree."""
        if not self.boundary.contains(point):
            # The point does not lie inside boundary: bail.
            return False

        if not self.divided:
          if self.point == None:
              # Node doesn't have a point yet.
              self.point = point
              return True

          # Already leaf: divide if necessary, then try the sub-quads.
          self.divide()

        return (self.ne.insert(point) or
                self.nw.insert(point) or
                self.se.insert(point) or
                self.sw.insert(point))

    def get_points_rec(self, found_points):
        """Find the points in the quadtree that lie within boundary."""
        if self.point != None:
            found_points.append(self.point)

        # if this node has children, search them too.
        if self.divided:
          self.nw.get_points_rec(found_points)
          self.ne.get_points_rec(found_points)
          self.se.get_points_rec(found_points)
          self.sw.get_points_rec(found_points)
        return found_points

    def get_points(self):
      return self.get_points_rec([])

    def __len__(self):
        """Return the number of points in the quadtree."""
        if self.point != None:
          npoints = 1
        else:
          npoints = 0
        if self.divided:
            npoints += len(self.nw)+len(self.ne)+len(self.se)+len(self.sw)
        return npoints

# POINT QUADTREE ##############################################################
class PointQuadTree:
    """Point Quadtree implementation."""

    def __init__(self, boundary, ax, depth=0):
        """Initialize this node of the quadtree."""
        self.boundary = boundary # boundary of current block
        self.ax = ax
        self.depth = depth # mostly for string visualization spacing
        self.point = None # center point
        self.connection = [] # WSP connection
        self.divided = False # flag for if divided into 4 child quads

    def __str__(self):
        """Return a string representation of this node, suitably formatted."""
        sp = ' ' * self.depth * 2
        s = str(self.boundary) + ' --> ' + str(self.point) 
        if not self.divided:
            return s
        return s + '\n' + '\n'.join([
                sp + 'nw: ' + str(self.nw), sp + 'ne: ' + str(self.ne),
                sp + 'se: ' + str(self.se), sp + 'sw: ' + str(self.sw)])

    def str_short(self):
      return str(self.boundary)

    def diameter(self):
      return self.boundary.diameter()
    def center(self):
      return self.boundary.center()

    def divide(self):
        """Divide (branch) this node by spawning four children nodes around a point."""
        self.nw = PointQuadTree(Rect(self.boundary.xMin, self.point.y, self.point.x, self.boundary.yMax), self.ax, self.depth + 1)
        self.ne = PointQuadTree(Rect(self.point.x, self.point.y, self.boundary.xMax, self.boundary.yMax), self.ax, self.depth + 1)
        self.se = PointQuadTree(Rect(self.point.x, self.boundary.yMin, self.boundary.xMax, self.point.y), self.ax, self.depth + 1)
        self.sw = PointQuadTree(Rect(self.boundary.xMin, self.boundary.yMin, self.point.x, self.point.y), self.ax, self.depth + 1)
        self.divided = True
        # draw
        self.ax[0].plot([self.point.x, self.point.x],[self.boundary.yMin, self.boundary.yMax], color="gray")
        self.ax[0].plot([self.boundary.xMin, self.boundary.xMax],[self.point.y, self.point.y], color="gray")
        self.ax[1].plot([self.point.x, self.point.x],[self.boundary.yMin, self.boundary.yMax], color="gray")
        self.ax[1].plot([self.boundary.xMin, self.boundary.xMax],[self.point.y, self.point.y], color="gray")

    def insert(self, point):
        """Try to insert Point point into this QuadTree."""
        if not self.boundary.contains(point):
            # The point does not lie inside boundary: bail.
            return False

        if self.point == None:
            # Node doesn't have a point yet.
            self.point = point
            return True

        # Already leaf: divide if necessary, then try the sub-quads.
        if not self.divided:
            self.divide()

        return (self.ne.insert(point) or
                self.nw.insert(point) or
                self.se.insert(point) or
                self.sw.insert(point))

    def get_points_rec(self, found_points):
        """Find the points in the quadtree that lie within boundary."""
        if self.point != None:
            found_points.append(self.point)

        # if this node has children, search them too.
        if self.divided:
          self.nw.get_points_rec(found_points)
          self.ne.get_points_rec(found_points)
          self.se.get_points_rec(found_points)
          self.sw.get_points_rec(found_points)
        return found_points

    def get_points(self):
      return self.get_points_rec([])

    def __len__(self):
        """Return the number of points in the quadtree."""
        if self.point != None:
          npoints = 1
        else:
          npoints = 0
        if self.divided:
            npoints += len(self.nw)+len(self.ne)+len(self.se)+len(self.sw)
        return npoints