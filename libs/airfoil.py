# coding=utf-8

import numpy as np
import sympy as sym
from math import *

### BLADE: BEGIN ###

class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        
class PrimitivePoints:
     def __init__(self, x_coordinates = np.array([]), y_coordinates = np.array([]), z_coordinates = np.array([])): #x_coordinates - np.array, #y_coordinates - np.array
        self.points = np.array([])
        for i in range(len(x_coordinates)):
            self.points = np.append(self.points, Point(x_coordinates[i], y_coordinates[i], z_coordinates[i]))
     def addPoint(self, point):
        self.points = np.append(self.points, point)
     def setPoints(self, points):
        self.points = np.array(points)
        
def createCirclePoints(centerPoint, R, pointsCount): #centerPoint - Point, #R - double, #pointsCount - int 
    t_array = np.linspace(0, 2*pi, pointsCount)    
    y_array = np.array([])
    x_array = np.array([])
    for t in t_array:
        x_array = np.append(x_array, centerPoint.x + R*cos(t))
        y_array = np.append(y_array, centerPoint.y + R*sin(t))
    z_array = np.full(len(t_array),centerPoint.z) 
    return PrimitivePoints(x_array, y_array, z_array)
    
def createEdgePoints(edgeCircle, pt1, pt2, edgeType):
    k = (pt2.y-pt1.y)/(pt2.x-pt1.x)
    if ((k > 0 and edgeType == 't') or (k < 0 and edgeType == 'l')):     
        edgeCircle.points = filter(lambda point: (point.y < (k*(point.x-pt1.x)+pt1.y)), edgeCircle.points)
    else:
        edgeCircle.points = filter(lambda point: (point.y > (k*(point.x-pt1.x)+pt1.y)), edgeCircle.points)
    edgeCircle.points = sorted(edgeCircle.points, key=lambda point: point.y)
    return edgeCircle
    
def drawPrimitivePoints(geo, primitivePoints):
    for i in range(len(primitivePoints.points)):
        point = geo.createPoint()
        point.setPosition((primitivePoints.points[i].x, primitivePoints.points[i].y, primitivePoints.points[i].z))
        
def drawPoint(geo, point):
    p = geo.createPoint()
    p.setPosition((point.x, point.y, point.z))
    
def drawRectangle(geo, p1, p2, p3, p4): 
    poly = geo.createPolygon()
    poly.addVertex(p1)
    poly.addVertex(p2)
    poly.addVertex(p3)
    poly.addVertex(p4)
        
def touchPointsLineCircle(linePoint, circleCenter, R):
    
    xc = circleCenter.x
    yc = circleCenter.y    
    xl = linePoint.x
    yl = linePoint.y
    h = linePoint.z
    
    k1 = (-(xl-xc)*(yl-yc)+R*sqrt((xc-xl)**2+(yc-yl)**2-R**2))/(R**2-(xc-xl)**2)    
    k2 = (-(xl-xc)*(yl-yc)-R*sqrt((xc-xl)**2+(yc-yl)**2-R**2))/(R**2-(xc-xl)**2) 
    a1 = k1**2+1
    a2 = k2**2+1
    b1 = -2*xc-2*k1*(yc-yl+k1*xl)
    b2 = -2*xc-2*k2*(yc-yl+k2*xl)   
    #c = (yc-yl+k*xl)**2-R**2+xc**2
        
    xt = -b1/(2*a1)
    yt = k1*(xt-xl)+yl
    firstTouchPoint = Point(xt, yt, h)
    
    xt = -b2/(2*a2)
    yt = k2*(xt-xl)+yl
    secondTouchPoint = Point(xt, yt, h)
    
    return np.array([firstTouchPoint, secondTouchPoint])

### BLADE: END ###

### FRENE FRAME: BEGIN ###

a0,a1,a2,x,y,z,t = sym.symbols('a0 a1 a2 x y z t')

def getTNBVectorsByVectorFuntion(x,y,z):
    alphaSym = [x, y, z]; dAlphaSym = []; ddAlphaSym = []
    for coorFunSym in alphaSym:
        dAlphaSym.append(sym.diff(coorFunSym, t))
    for dCoorFunSym in dAlphaSym:
        ddAlphaSym.append(sym.diff(dCoorFunSym, t)) 
    tauSym = []
    mdAlphaSym = sym.sqrt(dAlphaSym[0]**2+dAlphaSym[1]**2+dAlphaSym[2]**2)    
    for dCoorFunSym in dAlphaSym: tauSym.append(dCoorFunSym/mdAlphaSym)
    bnormSym = []
    bnormNominatorSym = np.cross(dAlphaSym,ddAlphaSym)
    mbnormSym = sym.sqrt(bnormNominatorSym[0]**2+bnormNominatorSym[1]**2+bnormNominatorSym[2]**2)
    for coorFunSym in bnormNominatorSym: bnormSym.append(coorFunSym/mbnormSym)
    normSym = np.cross(bnormSym,tauSym)
    normSym = normSym.tolist()
    return [tauSym, bnormSym]

def getPolynom2DCoefficients(x0, x1, x2, y0, y1, y2):
    A = np.array([[x0**2,x0,1],[x1**2,x1,1],[x2**2,x2,1]])
    B = np.array([y0,y1,y2])
    C = np.linalg.solve(A,B)   
    return C  
    
def getFreneFrame(u1,u2,v1,v2,p):
    xu = t;                 xv = v1.position().x() 
    yu = a0*t**2+a1*t+a2;   yv = a0*t**2+a1*t+a2
    zu = u1.position().z(); zv = t
    [tauSym, bnormSym] = getTNBVectorsByVectorFuntion(xu,yu,zu)
    Cu = getPolynom2DCoefficients(u1.position().x(), p.position().x(), u2.position().x(), u1.position().y(), p.position().y(), u2.position().y())
    Cv = getPolynom2DCoefficients(v1.position().z(), p.position().z(), v2.position().z(), v1.position().y(), p.position().y(), v2.position().y())
    tau = []; bnorm = []
    for coor_fun in tauSym: tau.append(float(coor_fun.subs([(t,p.position().x()),(a0,Cu[0]),(a1,Cu[1]),(a2,Cu[2])])))
    for coor_fun in bnormSym: bnorm.append(float(coor_fun.subs([(t,p.position().z()),(a0,Cv[0]),(a1,Cv[1]),(a2,Cv[2])])))
    norm = np.cross(tau,bnorm)
    tau = np.array(tau); tau = tau.reshape((-1,1))
    bnorm = np.array(bnorm); bnorm = bnorm.reshape((-1,1))
    norm = norm.reshape((-1,1))
    # we get tau -> x; bnorm -> y, norm -> z because norm(z) = tau(x) x bnorm(y) => so we need [tau, bnorm, norm] matrix instead of [tau, norm, bnorm]
    return np.hstack((tau, np.hstack((bnorm,norm))))  
    
def drawFreneFrame(geo, p, tau, bnorm, norm):
    tauPoint = geo.createPoint()
    tauPoint.setPosition((float(p.position().x()+tau[0]), float(p.position().y()+tau[1]), float(p.position().z()+tau[2])))
    tauPoly = geo.createPolygon()
    tauPoly.addVertex(p)
    tauPoly.addVertex(tauPoint)
    
    bnormPoint = geo.createPoint()
    bnormPoint.setPosition((float(p.position().x()+bnorm[0]), float(p.position().y()+bnorm[1]), float(p.position().z()+bnorm[2])))
    normPoly = geo.createPolygon()
    normPoly.addVertex(p)
    normPoly.addVertex(bnormPoint)
    
    normPoint = geo.createPoint()
    normPoint.setPosition((float(p.position().x()+norm[0]), float(p.position().y()+norm[1]), float(p.position().z()+norm[2])))
    normpPoly = geo.createPolygon()
    normpPoly.addVertex(p)
    normpPoly.addVertex(normPoint)
    
### FRENE FRAME: END ###