# coding=utf-8

import numpy as np
from abc import ABCMeta, abstractmethod
from enum import Enum

def shift(arr, num, fill_value):
    result = np.empty_like(arr)
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result

class CurveType(Enum):
    LIN = 1
    CIRC = 2
    CUBIC_SPLINE = 3
    QUADRATIC_SPLINE = 4
    POLYNOMIAL = 5

 # max distance between trajectory interior points 1mm, because max velocity 250 mm/s
 # and time limit for data frame exchange is 4ms => 4*10^-3*250 = 1mm  
MAX_DIS = 1

### ------------------------------ CURVE ------------------------------ ###
class InterpolationCurve:
    def __init__(self, nodePoints, interpolationStrategy, step):
        self.__nodePoints = nodePoints 
        self.__interpolationStrategy = interpolationStrategy
        self.__type = type
        self.__step = step

    @property
    def nodePoints(self):
        return self.__nodePoints
    @property
    def interpolationStrategy(self):
        return self.__interpolationStrategy
    @property
    def step(self):
        return self.__step

    @nodePoints.setter
    def setNodePoints(self,nodePoints):
        self.__nodePoints = nodePoints
    @interpolationStrategy.setter
    def setInterpolationStrategy(self,interpolationStrategy):
        self.__interpolationStrategy = interpolationStrategy
    @step.setter
    def setStep(self,step):
        self.__step = step

    def interpolate(self):
        self.__interpolationStrategy.execute(self.__nodePoints, self.__step)
          
 #   def drawCurve(): # just create and connect all the points using addVertex  
 #       print("draw")
### ------------------------------ CURVE ------------------------------ ###

### ------------------------------------------------------------------- ###
# Strategy for interpolation
class InterpolationStrategy: #return list of points along the curve
    __metaclass__ = ABCMeta
    @abstractmethod
    def execute(self, nodePoints, step):
        pass
### ------------------------------------------------------------------- ###

class LinearInterpolation(InterpolationStrategy):
    def execute(self, nodePoints, step):
        n = len(nodePoints)
        h = np.empty([n-1,2,3])
        for i in range(n-1):
            h[i] = np.array([nodePoints[i+1]-nodePoints[i],nodePoints[i]])
            
        # draw polygon
        t_array = np.arange(0.0, 1.0, step)
        for i in range(n-1):
            for t in t_array:
                p = geo.createPoint()
                pos = h[i,0]*t + h[i,1]
                p.setPosition((pos[0],pos[1],pos[2]))

class CircInterpolation(InterpolationStrategy):
    def execute(self, nodePoints, step):
        print("circ")

class QuadraticInterpolation(InterpolationStrategy):
    def execute(self, nodePoints, step):
        print("quad")

class CubicInterpolation(InterpolationStrategy):
    def execute(self, nodePoints, step):
        n = len(nodePoints)
        M = np.empty([n-2,n-2])
        M[0] = np.concatenate((np.array([7,2]),np.zeros(n-4)));
        M[-1] = np.concatenate((np.zeros(n-4),np.array([2,7])))

        P = np.empty([n-2,3])
        P[0] = -3*(nodePoints[0]+nodePoints[1]-2*nodePoints[2])
        P[-1] = -3*(2*nodePoints[-3]-nodePoints[-2]-nodePoints[-1])

        row = np.concatenate((np.array([1,4,1]),np.zeros(n-2)))

        for i in range(1,n-3):
            M[i] = row
            P[i] = np.array([3*(nodePoints[i+2]-nodePoints[i])])
            row = shift(row, 1, 0)
            
        dP = np.linalg.solve(M,P) # dP[0] = dP2, dP[1] = dP3, dP[2] = dP4 ...

        # Calculate dP0 and dPn
        dP0 = np.array([1.5*(nodePoints[1]-nodePoints[0])-0.5*dP[0]])
        dPn = np.array([1.5*(nodePoints[-1]-nodePoints[-2])-0.5*dP[-1]])
        dP = np.concatenate((dP0,dP), axis = 0)
        dP = np.concatenate((dP,dPn), axis = 0)

        # Hermit Matrix
        H = np.array([[2,-2,1,1],[-3,3,-2,-1],[0,0,1,0],[1,0,0,0]])

        # Calculate Hermite Splines cofficients h[i], starting from 0...n-1
        # every hermite spline represented by 4 coefficients 
        # every coefficient is a row-vector 1x3, so we need matrix with dim = (n-1)x4x3
        h = np.empty([n-1,4,3])
        for i in range(n-1):
            h[i] = np.dot(H,np.array([nodePoints[i],nodePoints[i+1],dP[i],dP[i+1]]))
            
        # draw cubic spline
        t_array = np.arange(0.0, 1.0, step)
        for i in range(n-1):
            for t in t_array:
                p = geo.createPoint()
                pos = h[i,0]*t**3 + h[i,1]*t**2 + h[i,2]*t + h[i,3]
                p.setPosition((pos[0],pos[1],pos[2]))

class PolyInterpolation(InterpolationStrategy):
    def execute(self, nodePoints, step):
        print("poly")
        




   
