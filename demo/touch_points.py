import numpy as np
from dataclasses import dataclass
import math

### A4-A4 ###

@dataclass
class Point:
    x: float
    y: float

@dataclass
class Circle:
    center: Point
    radius: float

def find_tangents(circle: Circle, p: Point):
    dx, dy = p.x - circle.center.x, p.y - circle.center.y
    A = dx**2 + dy**2
    s = math.sqrt(A - circle.radius**2)
    denom = A
    t1 = Point(circle.center.x + (circle.radius**2 * dx + circle.radius * dy * s) / denom,
               circle.center.y + (circle.radius**2 * dy - circle.radius * dx * s) / denom)
    t2 = Point(circle.center.x + (circle.radius**2 * dx - circle.radius * dy * s) / denom,
               circle.center.y + (circle.radius**2 * dy + circle.radius * dx * s) / denom)
    return t1, t2

def all_tangents(circle: Circle, p1: Point, p2: Point):
    return find_tangents(circle, p1), find_tangents(circle, p2)

# Входная кромка R1

p1_cxe = [-11.5, 0.45]
p1_cve = [-11.5, -0.08]
p1_c = [-11.76, 0.1]
R1 = 0.2

# Выходная кромка R2

p2_cxe = [12, 0.65]
p2_cve = [12, 0.18]
p2_c = [13.51, 0.17]
R2 = 0.14

# -----

circle1 = Circle(Point(p1_c[0], p1_c[1]), R1)
p1, p2 = Point(p1_cxe[0], p1_cxe[1]), Point(p1_cve[0], p1_cve[1])
tangents = all_tangents(circle1, p1, p2)

for tangent in tangents:
    for t in tangent:
        print(t)

# p1_cxt = (-11.847953, 0.279622)
# p1_cvt = (-11.744182, -0.099373)

# p2_cxt = (13.540457, 0.306647)
# p2_cvt = (13.496097, 0.030692)

