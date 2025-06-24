import numpy as np
import json
from linalg import *
from diffgeom import *

def get_frene(ptc: np.array, pt: np.array):
    norm = (pt - ptc)/np.linalg.norm(pt - ptc) # normal
    tau = [-norm[1], norm[0], norm[2]]
    binorm = np.cross(tau, norm)
    return Frene(tau, norm, binorm, pt)

with open("99.01.25.242.json", 'r') as file:
    blade = json.load(file)

pt1 = np.array([-11.655, -0.172, 61.8])
pt2 = np.array([-11.683, -0.317, 61.8])
pt3 = np.array([-11.5, -0.42, 61.8])

print(points2circ(pt1, pt2, pt3))

