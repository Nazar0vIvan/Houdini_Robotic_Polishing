import numpy as np
from math import *
from linalg import *

def plane(norm: list, pt: list) -> list:
    A = norm[0]
    B = norm[1]
    C = norm[2]
    D = -A*pt[0]-B*pt[1]-C*pt[2]
    return np.array([A,B,C,D])

# {2} - scheme
# {1} - root
# {0} - airfoil
# {2} == {1}

np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)

pt1_2 = np.array([[-9.9], [8.44], [9.6]])
pt2_2 = np.array([[12.2], [8.44], [14.3]])

k = (pt2_2[2][0] - pt1_2[2][0])/(pt2_2[0][0] - pt1_2[0][0])
b = pt1_2[2][0] - k*pt1_2[0][0]

norm_2 = normalize(-np.array([[k],[0.],[-1]]))

#print(k)
#print(b)
#print(norm_2)

shelf_plane = plane(norm_2.flatten(), pt2_2.flatten())
print(shelf_plane)
'''
R20 =  rotationMatrix3x3(radians(-41.0), 'z') # {2} -> {0}

pt1_0 = R20 @ pt1_2
pt2_0 = R20 @ pt2_2

[k,b] = poly((pt1_0[0][0],pt1_0[2][0]),(pt2_0[0][0],pt2_0[2][0]))

norm_0 = normalize(-np.array([[k],[0.],[-1]]))
shelf_plane = plane(norm_0.flatten(), pt2_0.flatten())
print(shelf_plane)

R02 = rotationMatrix3x3(radians(41.0), 'z')
print(R02)
'''
# check