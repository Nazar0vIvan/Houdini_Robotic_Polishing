import numpy as np
import sympy as sym
import hou.session as hs 
#import sys, os
#sys.path.append(os.environ["HIP"]+"/libs")

from airfoil import *

node = hou.pwd()
geo = node.geometry()

# leading edge - входная кромка (тольще)
R2 = np.array([ 0.10,  0.10,  0.10,  0.12,  0.14,  0.16,  0.19,  0.21,  0.23,  0.25])
X2 = np.array([13.32, 13.47, 13.48, 13.51, 13.51, 13.50, 13.45, 13.41, 13.37, 13.34])
Y2 = np.array([ 3.14,  2.10,  1.88,  1.12,  0.17, -0.72, -1.55, -1.96, -2.34, -2.70])

# trailing edge - выходная кромка (тоньше)
R1 = np.array([  0.09,   0.11,   0.12,   0.15,   0.20,   0.26,   0.33,   0.38,   0.43,   0.49])
X1 = np.array([-11.79, -11.80, -11.80, -11.80, -11.76, -11.69, -11.60, -11.55, -11.49, -11.43])
Y1 = np.array([ -0.78,  -0.37,  -0.35,  -0.08,   0.11,   0.14,   0.08,   0.02,  -0.07,  -0.19])

X = np.concatenate((np.arange(-11.5,-7.75,0.25), np.arange(-6,14,2)))
# Use drop down menu to select examples.

#concave - корыто
Y_cv = np.array([[-0.78, -0.70, -0.62, -0.54, -0.45, -0.37, -0.29, -0.20, -0.12, -0.04,  0.04,  0.12,  0.20,  0.28,  0.35,  0.90,  1.38,  1.78,  2.11,  2.38,  2.60,  2.77,  2.90,  2.98,  3.03], #A1-A1
                 [-0.42, -0.36, -0.30, -0.24, -0.18, -0.12, -0.06,  0.01,  0.07,  0.13,  0.19,  0.25,  0.31,  0.37,  0.43,  0.84,  1.19,  1.47,  1.69,  1.86,  1.98,  2.06,  2.09,  2.09,  2.05], #A2-A2
                 [-0.40, -0.35, -0.29, -0.24, -0.18, -0.12, -0.06,  0.00,  0.06,  0.12,  0.17,  0.23,  0.28,  0.34,  0.40,  0.78,  1.11,  1.37,  1.57,  1.72,  1.82,  1.88,  1.90,  1.89,  1.83], #A2a-A2a
                 [-0.19, -0.16, -0.12, -0.08, -0.04,  0.00,  0.05,  0.09,  0.14,  0.18,  0.22,  0.27,  0.31,  0.35,  0.39,  0.68,  0.92,  1.09,  1.21,  1.29,  1.33,  1.33,  1.29,  1.21,  1.10], #A3-A3
                 [-0.08, -0.06, -0.05, -0.03,  0.00,  0.03,  0.06,  0.08,  0.11,  0.14,  0.17,  0.20,  0.23,  0.26,  0.29,  0.48,  0.62,  0.71,  0.75,  0.74,  0.70,  0.63,  0.51,  0.37,  0.18], #A4-A4
                 [-0.12, -0.12, -0.12, -0.11, -0.10, -0.08, -0.06, -0.05, -0.03, -0.01,  0.01,  0.03,  0.05,  0.07,  0.09,  0.21,  0.28,  0.30,  0.27,  0.21,  0.10, -0.04, -0.21, -0.42, -0.67], #A5-A5
                 [-0.26, -0.28, -0.29, -0.30, -0.30, -0.29, -0.28, -0.27, -0.26, -0.25, -0.24, -0.23, -0.22, -0.21, -0.19, -0.13, -0.11, -0.14, -0.23, -0.35, -0.51, -0.70, -0.93, -1.20, -1.50], #A6-A6
                 [-0.39, -0.41, -0.42, -0.43, -0.43, -0.43, -0.42, -0.42, -0.41, -0.40, -0.39, -0.39, -0.38, -0.37, -0.33, -0.33, -0.39, -0.50, -0.64, -0.82,  1.04, -1.29, -1.58, -1.91]])       #A7-A7

#[0.24,  0.36,  0.47,  0.57, 0.66, 0.74, 0.83, 0.92, 1.00, 1.07, 1.14, 1.21, 1.28, 1.35, 1.41, 1.83, 2.12, 2.31, 2.41, 2.43, 2.36, !2.23!,  !2.03!,  1.79,  1.49], #A3-A3
#[ 0.52,  0.65,  0.76,  0.86, 0.95, 1.03, 1.11, 1.19, 1.27, 1.33, 1.39, 1.45, 1.51, 1.57, 1.62, 1.93, 2.11, 2.16, 2.11, 1.95, 1.69, 1.35,  0.93,  0.44, -0.11], #A5-A5
#convex - спинка
Y_cx = np.array([[-0.50, -0.36, -0.23, -0.10, 0.02, 0.13, 0.24, 0.35, 0.46, 0.56, 0.66, 0.76, 0.85, 0.95, 1.03, 1.67, 2.17, 2.58, 2.89, 3.12, 3.27, 3.36,  3.39,  3.37,  3.30], #A1-A1
                 [-0.09,  0.04,  0.15,  0.27, 0.36, 0.46, 0.56, 0.65, 0.75, 0.83, 0.91, 0.99, 1.07, 1.15, 1.22, 1.73, 2.12, 2.42, 2.61, 2.73, 2.77, 2.75,  2.67,  2.54,  2.36], #A2-A2
                 [-0.05,  0.08,  0.19,  0.30, 0.39, 0.49, 0.58, 0.68, 0.77, 0.85, 0.92, 1.00, 1.08, 1.16, 1.22, 1.72, 2.09, 2.36, 2.53, 2.63, 2.65, 2.61,  2.51,  2.35,  2.15], #A2a-A2a
                 [ 0.24,  0.36,  0.47,  0.57, 0.66, 0.74, 0.83, 0.92, 1.00, 1.07, 1.14, 1.21, 1.28, 1.35, 1.41, 1.83, 2.12, 2.31, 2.41, 2.43, 2.36, 2.23,  2.03,  1.79,  1.49], #A3-A3
                 [ 0.45,  0.57,  0.68,  0.78, 0.86, 0.94, 1.02, 1.10, 1.18, 1.24, 1.30, 1.37, 1.43, 1.49, 1.54, 1.89, 2.11, 2.23, 2.24, 2.16, 1.99, 1.76,  1.44,  1.08,  0.65], #A4-A4
                 [ 0.52,  0.65,  0.76,  0.86, 0.95, 1.03, 1.11, 1.19, 1.27, 1.33, 1.39, 1.45, 1.51, 1.57, 1.62, 1.93, 2.11, 2.16, 2.11, 1.95, 1.69, 1.35,  0.93,  0.44, -0.11], #A5-A5
                 [ 0.52,  0.65,  0.78,  0.89, 0.98, 1.06, 1.15, 1.23, 1.31, 1.38, 1.44, 1.50, 1.56, 1.62, 1.67, 1.98, 2.13, 2.14, 2.03, 1.80, 1.46, 1.01,  0.48, -0.13, -0.82], #A6-A6
                 [ 5.64,  0.77,  0.90,  0.99, 1.07, 1.16, 1.25, 1.34, 1.41, 1.47, 1.53, 1.59, 1.66, 1.71, 2.02, 2.16, 2.16, 2.02, 1.75, 1.37, 0.87, 0.28, -0.40, -1.16]])       #A7-A7
     
Z = np.array([np.full(len(X),71.8), 
              np.full(len(X),61.8),
              np.full(len(X),59.8),
              np.full(len(X),51.8),
              np.full(len(X),41.8),
              np.full(len(X),31.8),
              np.full(len(X),21.8),
              np.full(len(X),16.8),
              np.full(len(X),11.8),
              np.full(len(X),6.8)])     
                 
profilesCount = len(Y_cv)

# CROSS SECTIONS
for i in range(1, profilesCount-1): # , A2-A2 to A6-A6 
    concave = PrimitivePoints(X, Y_cv[i], Z[i])
    convex = PrimitivePoints(X, Y_cx[i], Z[i])
    convex.points = np.flip(convex.points, 0)
        
    # draw CONCAVE points
    n1 = len(geo.points())
    drawPrimitivePoints(geo, concave)
    k1 = len(geo.points())
    #curConcavePointGroupName = "A"+str(i)+"cv"
    #curConcavePointGroup = geo.createPointGroup(curConcavePointGroupName)
    #curConcavePointGroup.add(geo.points()[n:k])
    # draw CONCAVE polygon
    profile = geo.createPolygon()
    profile.setIsClosed(False)
    for point in geo.points()[n1:k1]:
        profile.addVertex(point)
        
    # draw CONVEX points
    n2 = len(geo.points())    
    drawPrimitivePoints(geo, convex)
    k2 = len(geo.points())
    #curConcavePointGroupName = "A"+str(i)+"cx"
    #curConcavePointGroup = geo.createPointGroup(curConcavePointGroupName)
    #curConcavePointGroup.add(geo.points()[n:k])
    # draw CONVEX polygon
    profile = geo.createPolygon()
    profile.setIsClosed(False)    
    for point in geo.points()[n2:k2]:
        profile.addVertex(point)

'''
    leadingEdge = geo.createPolygon()
    leadingEdge.setIsClosed(True)
    leadingEdge.addVertex(geo.points()[k2-1])
    leadingEdge.addVertex(geo.points()[n1])
    
    trailingEdge = geo.createPolygon()
    trailingEdge.setIsClosed(True)
    trailingEdge.addVertex(geo.points()[k1-1])
    trailingEdge.addVertex(geo.points()[n2])
'''
# AIRFOIL SURFACE




        