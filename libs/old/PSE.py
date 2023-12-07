import hou
import numpy as np
import sympy as sym
from math import *
import sys

from abc import ABCMeta, abstractmethod
from enum import Enum

deg = 180/pi

def rotationMatrix(angle, axis):
    if (axis == "x"):
        rot = np.array([[1.0,0.0,0.0],[0.0,np.cos(angle),-np.sin(angle)],[0.0,sin(angle),np.cos(angle)]])        
    elif (axis == "y"):
        rot = np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-sin(angle),0,cos(angle)]])
    elif (axis == "z"):
        rot = np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])       
    rot[np.absolute(rot)<=0.0001] = 0.0
    return rot

def translationMatrix(dx,dy,dz):
    return np.array([[1,0,0,dx],[0,1,0,dy],[0,0,1,dz],[0,0,0,1]])

def eulerAnglesToRotationMatrix(A = 0, B = 0, C = 0):
    A = A*pi/180.0; B = B*pi/180.0; C = C*pi/180.0
    result = np.dot(np.dot(rotationMatrix(A,"z"),rotationMatrix(B,"y")),rotationMatrix(C,"x"))
    result[np.absolute(result)<=0.0001]=0
    return result

def rotationMatrixToEulerAngles(rot):
    B1 = -np.arcsin(rot[2][0])
    B2 = pi + np.arcsin(rot[2][0])
    C1 = np.arctan2(rot[2][1]/np.cos(B1), rot[2][2]/np.cos(B1))
    C2 = np.arctan2(rot[2][1]/np.cos(B2), rot[2][2]/np.cos(B2))
    A1 = np.arctan2(rot[1][0]/np.cos(B1), rot[0][0]/np.cos(B1))
    A2 = np.arctan2(rot[1][0]/np.cos(B2), rot[0][0]/np.cos(B2))  
    return {'A1':A1, 'A2':A2, 'B1':B1, 'B2':B2, 'C1':C1, 'C2':C2}

def dhMatrix(q,i):
    qf = q + q0[i];
    mat = np.array([[cos(qf),               -sin(qf),                0,              a[i]],
                    [cos(alpha[i])*sin(qf),  cos(alpha[i])*cos(qf), -sin(alpha[i]), -d[i]*sin(alpha[i])],
                    [sin(alpha[i])*sin(qf),  sin(alpha[i])*cos(qf),  cos(alpha[i]),  d[i]*cos(alpha[i])],
                    [0,                      0,                      0,              1]]);
    return mat

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

def points2plane(x,y,z):
    n = len(x);

    U = np.zeros((3,3))
    V = np.zeros(3)
    
    U[0,0] = np.sum(x**2)
    U[0,1] = np.sum(x*y)
    U[0,2] = np.sum(x)
    
    U[1,0] = U[0,1];
    U[1,1] = np.sum(y**2)
    U[1,2] = np.sum(y);
    
    U[2,0] = U[0,2];
    U[2,1] = U[1,2];
    U[2,2] = n;
    
    V[0] = np.sum(x*z)
    V[1] = np.sum(y*z)
    V[2] = sum(z)
    
    P = np.linalg.solve(U,V)

    # AA*x + BB*y - z + DD = 0
    AA = P[0]
    BB = P[1]
    DD = P[2]
    
    # A*x + B*y + C*z + D = 0 => z = -(A/C)*x - (B/C)*y - D/C
    # z = AA*x + BB*y + DD => AA = -(A/C); BB = -(B/C); DD = -D/C
    
    # A^2+B^2+C^2 = 1
    # A = -AA*C; -B = BB*C; D = -DD*C
    # (-AA*C)^2 + (-BB*C)^2 + C^2 = 1 =>
    C = np.sqrt(1/(AA**2+BB**2+1))
    A = -AA*C;
    B = -BB*C;
    D = -DD*C;
    
    return [A,B,C,D,AA,BB,DD]

### FRENE FRAME: BEGIN ###

a0,a1,a2,x,y,z,t = sym.symbols('a0 a1 a2 x y z t')

def getTNBVectorsByVectorFunction(x,y,z):
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
    [tauSym, bnormSym] = getTNBVectorsByVectorFunction(xu,yu,zu)
    Cu = getPolynom2DCoefficients(u1.position().x(), p.position().x(), u2.position().x(), u1.position().y(), p.position().y(), u2.position().y())
    Cv = getPolynom2DCoefficients(v1.position().z(), p.position().z(), v2.position().z(), v1.position().y(), p.position().y(), v2.position().y())
    tau = []; bnorm = []s
    for coor_fun in tauSym: tau.append(float(coor_fun.subs([(t,p.position().x()),(a0,Cu[0]),(a1,Cu[1]),(a2,Cu[2])])))
    for coor_fun in bnormSym: bnorm.append(float(coor_fun.subs([(t,p.position().z()),(a0,Cv[0]),(a1,Cv[1]),(a2,Cv[2])])))
    norm = np.cross(tau,bnorm)
    tau = np.array(tau); tau = tau.reshape((-1,1))
    bnorm = np.array(bnorm); bnorm = bnorm.reshape((-1,1))
    norm = norm.reshape((-1,1))
    # we get tau -> x; bnorm -> y, norm -> z because norm(z) = tau(x) x bnorm(y) => so we need [tau, bnorm, norm] matrix instead of [tau, norm, bnorm]
    freneFrame = np.hstack((tau, np.hstack((bnorm,norm))))
    return freneFrame if norm[1] > 0 else np.dot(freneFrame,np.array([[1,0,0],[0,-1,0],[0,0,-1]]))

def drawFreneFrame(geo, p, tau, bnorm, norm, scale=10):
    tauPoint = geo.createPoint()
    tauPoint.setPosition((float(p.position().x()+scale*tau[0]), float(p.position().y()+scale*tau[1]), float(p.position().z()+scale*tau[2])))
    tauPoly = geo.createPolygon()
    tauPoly.addVertex(p)
    tauPoly.addVertex(tauPoint)
    
    bnormPoint = geo.createPoint()
    bnormPoint.setPosition((float(p.position().x()+scale*bnorm[0]), float(p.position().y()+scale*bnorm[1]), float(p.position().z()+scale*bnorm[2])))
    bnormPoly = geo.createPolygon()
    bnormPoly.addVertex(p)
    bnormPoly.addVertex(bnormPoint)
    
    normPoint = geo.createPoint()
    normPoint.setPosition((float(p.position().x()+scale*norm[0]), float(p.position().y()+scale*norm[1]), float(p.position().z()+scale*norm[2])))
    normPoly = geo.createPolygon()
    normPoly.addVertex(p)
    normPoly.addVertex(normPoint)

### FRENE FRAME: END ###

