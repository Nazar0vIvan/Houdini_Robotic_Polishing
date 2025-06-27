import numpy as np
from dataclasses import dataclass

@dataclass
class Circle:
    center: np.ndarray
    radius: float

def find_furthest_tangents(p1: np.ndarray, p2: np.ndarray, circle: Circle) -> tuple[np.ndarray, np.ndarray]:
    pc, rc = circle.center, circle.radius

    def get_tangent(p, other_p):
        v = pc - p
        d_sq = np.dot(v, v)
        a = np.sqrt(max(d_sq - rc**2, 0))
        perp = np.array([-v[1], v[0]])
        t1 = p + (a**2 * v + a * rc * perp) / d_sq
        t2 = p + (a**2 * v - a * rc * perp) / d_sq
        return t1 if np.linalg.norm(t1 - other_p) > np.linalg.norm(t2 - other_p) else t2

    return get_tangent(p1, p2), get_tangent(p2, p1)

