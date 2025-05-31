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

np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)

pt1_2 = np.array([[-22.5], [-22.5*np.tan(radians(25.0))], [8.6]])
pt2_2 = np.array([[25.5],  [25.5*np.tan(radians(25.0))], [10.5]])

R21 = rotationMatrix3x3(radians(-25.0), 'z'); # {2} -> {1}

pt1_1 = R21 @ pt1_2
pt2_1 = R21 @ pt2_2

print(pt1_1)


[k,b] = poly((pt1_1[0][0],pt1_1[2][0]),(pt2_1[0][0],pt2_1[2][0]))

# norm_1 - нормаль к плоскости полки в СК root
norm_1 = normalize(-np.array([[k],[0.],[-1]]))

R10 = rotationMatrix3x3(radians(-17.0), 'z')

# norm_0 - нормаль к плоскости полки в СК airfoil
norm_0 = R10 @ norm_1
pt1_0 = R10 @ pt1_1
pt2_0 = R10 @ pt2_1

shelf_plane = plane(norm_0.flatten(), pt2_0.flatten())
#print(shelf_plane)

mat1 = np.array([[1,2,3],[4,5,6],[7,8,9]])
mat2 = np.array([[9,8,7],[6,5,4],[3,2,1]])

#print(mat1 @ mat2)

'''
R - матрица поворта из системы координат пера в систему координат, 
в которой Z направлена вдоль нормали к полке, а X и Y лежат в плоскости полки

norm = R * norm_0 -> [0,0,1]
'''
R = np.dot(rotationMatrix3x3(radians(17.0), 'z'), rotationMatrix3x3(np.arctan(k), 'y'))
#print(rotationMatrix3x3(np.arctan(k), 'y')) 


