import numpy as np

# check if edges intersect
def EdgesIntersect(e1, e2):
    a = e1.p1
    b = e1.p2
    c = e2.p1
    d = e2.p2
    cmp = Point(c.x - a.x, c.y - a.y)
    r = Point(b.x - a.x, b.y - a.y)
    s = Point(d.x - c.x, d.y - c.y)
    cmpxr = cmp.x * r.y - cmp.y * r.x
    cmpxs = cmp.x * s.y - cmp.y * s.x
    rxs = r.x * s.y - r.y * s.x
    if cmpxr == 0:
        return (c.x - a.x < 0) != (c.x - b.x < 0)
    if rxs == 0:
        return False
    rxsr = 1 / rxs
    t = cmpxs * rxsr
    u = cmpxr * rxsr
    return t >= 0 and t <= 1 and u >= 0 and u <= 1

# just the class for a 2D point
class Point:
    x = None
    y = None
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.arr = np.array([x,y])

# same as above, just the class for a 2D segment
class Segment:
    p1 = None
    p2 = None
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

# this is more important, a polygon can be stored that way.
class Polygon:
    points = None
    def __init__(self):
        self.points = []
    def AddPoint(self, p):
        self.points.append(p)
    def GetEdges(self):
        edges = []
        for i in range(len(self.points)):
            if i == len(self.points) - 1:
                i2 = 0
            else:
                i2 = i + 1
            edges.append(Segment(self.points[i], self.points[i2]))
        return edges
    def arraypoints(self):
        pts = self.points
        arraystuff = []
        for p in pts:
            arraystuff.append(p.arr)
        return np.asarray(arraystuff)
