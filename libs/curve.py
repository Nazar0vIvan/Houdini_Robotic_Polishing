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

### ------------------------------ CURVE ------------------------------ ###
# Strategy for interpolation
class InterpolationStrategy:
    __metaclass__ = ABCMeta
    def __init__(self):
        self.point = lambda i,t : None
    @abstractmethod
    def execute(self, nodePoints):
        pass
class LinearInterpolation(InterpolationStrategy):
    def execute(self, nodePoints):
        n = len(nodePoints)
        self.__coeffs = np.empty([n-1,2,3])
        for i in range(n-1):
            self.__coeffs[i] = np.array([nodePoints[i+1]-nodePoints[i],nodePoints[i]])
        self.point = lambda i,t : self.__coeffs[i,0]*t + self.__coeffs[i,1]
            
class CircInterpolation(InterpolationStrategy):
    def execute(self, nodePoints):
        print("circ interpolation done")

class QuadraticInterpolation(InterpolationStrategy):
    def execute(self, nodePoints):
        print("quadratic interpolation done")

class CubicInterpolation(InterpolationStrategy):
    def execute(self, nodePoints):
        n = len(nodePoints)
        M = np.empty([n-2,n-2])
        M[0] = np.concatenate((np.array([7,2]),np.zeros(n-4)));
        M[-1] = np.concatenate((np.zeros(n-4),np.array([2,7])))
        row = np.concatenate((np.array([1,4,1]),np.zeros(n-2)))

        P = np.empty([n-2,3])
        P[0] = -3*(nodePoints[0]+nodePoints[1]-2*nodePoints[2])
        P[-1] = -3*(2*nodePoints[-3]-nodePoints[-2]-nodePoints[-1])

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

        # Calculate Hermite Splines coefficients - coefficients[i], starting from 0...n-1
        # every hermite spline represented by 4 coefficients 
        # every coefficient is a row-vector 1x3, so we need matrix with dim = (n-1)x4x3
        self.__coeffs = np.empty([n-1,4,3])
        for i in range(n-1):
            P1 = nodePoints[i]
            P2 = nodePoints[i+1]
            dP1 = dP[i]
            dP2 = dP[i+1]
            self.__coeffs[i] = np.dot(H,np.array([P1,P2,dP1,dP2]))
        self.point = lambda i,t : self.__coeffs[i,0]*t**3 + self.__coeffs[i,1]*t**2 + self.__coeffs[i,2]*t + self.__coeffs[i,3]
#        print("cubic spline interpolation done")
            
class PolyInterpolation(InterpolationStrategy):
    def execute(self, nodePoints):
        print("polynomial interpolation done")

class InterpolationCurve:
    def __init__(self, nodePoints = np.array([]), interpolationStrategy = LinearInterpolation(), step = 0.2):
        self.__nodePoints = nodePoints 
        self.__interpolationStrategy = interpolationStrategy
        self.__step = step

    def nodePoints(self):
        return self.__nodePoints
    def points(self):
        return self.__points
    def interpolationStrategy(self):
        return self.__interpolationStrategy
    def step(self):
        return self.__step
    def setNodePoints(self,nodePoints):
        self.__nodePoints = nodePoints  
    def setInterpolationStrategy(self,interpolationStrategy):
        self.__interpolationStrategy = interpolationStrategy
    def setStep(self,step):
        self.__step = step
        
    def interpolate(self):
        if(self.__nodePoints.size != 0):
            self.__interpolationStrategy.execute(self.__nodePoints)
        else:
            hou.ui.displayMessage("There are no node points")
            
    def point(self,i,t):
        return self.__interpolationStrategy.point(i,t)
                  
    def calculateAndGetPoints(self): # just calculate all the curve points depending on self.__step value
        n = len(self.__nodePoints)
        self.__points = []
        t_array = np.arange(0.0, 1.0, self.__step)
        for i in range(n-1):
            for t in t_array:
                p = self.point(i,t)
                self.__points.append(np.array([p[0],p[1],p[2]]))
        self.__points = np.array(self.__points)
        return self.__points
### ------------------------------ CURVE ------------------------------ ###

        




   
