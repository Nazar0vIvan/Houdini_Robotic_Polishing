from dataclasses import dataclass
import sympy as sym
import numpy as np
from math import pi, sin, cos

def rotationMatrix3x3(angle: float, axis: float):
    if (axis == "x"):
        rot = np.array([[1.,0.,0.],[0.,np.cos(angle),-np.sin(angle)],[0.,sin(angle),np.cos(angle)]])        
    elif (axis == "y"):
        rot = np.array([[np.cos(angle),0.,np.sin(angle)],[0.,1.,0.],[-sin(angle),0.,cos(angle)]])
    elif (axis == "z"):
        rot = np.array([[np.cos(angle),-np.sin(angle),0.],[np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])       
    rot[np.absolute(rot)<=0.0001] = 0.
    return rot
    
def rotationMatrix4x4(angle: float, axis: float):
    if (axis == "x"):
        rot = np.array([[1.,0.,0.,0.],[0.,np.cos(angle),-np.sin(angle), 0.],[0.,sin(angle),np.cos(angle),0.],[0.,0.,0.,1.]])        
    elif (axis == "y"):
        rot = np.array([[np.cos(angle),0.,np.sin(angle),0.],[0.,1.,0.,0.],[-sin(angle),0.,cos(angle),0.],[0.,0.,0.,1.]])
    elif (axis == "z"):
        rot = np.array([[np.cos(angle),-np.sin(angle),0.,0.],[np.sin(angle),np.cos(angle),0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])       
    rot[np.absolute(rot)<=0.0001] = 0.0
    return rot

def translationMatrix(dx: float, dy: float, dz: float):
    return np.array([[1.,0.,0.,dx],[0.,1.,0.,dy],[0.,0.,1.,dz],[0.,0.,0.,1.]])

def euler2rot(A: float, B: float, C: float):
    A = A*pi/180.0; B = B*pi/180.0; C = C*pi/180.0
    result = np.dot(np.dot(rotationMatrix3x3(A, "z"), rotationMatrix3x3(B, "y")), rotationMatrix3x3(C, "x"))
    result[np.absolute(result)<=0.0001]=0
    return result

def rot2euler(rot, is_deg: bool = False):
    deg = 180/pi if is_deg else 1
    B1 = -np.arcsin(rot[2][0])
    B2 = pi + np.arcsin(rot[2][0])
    C1 = np.arctan2(rot[2][1]/np.cos(B1), rot[2][2]/np.cos(B1))
    C2 = np.arctan2(rot[2][1]/np.cos(B2), rot[2][2]/np.cos(B2))
    A1 = np.arctan2(rot[1][0]/np.cos(B1), rot[0][0]/np.cos(B1))
    A2 = np.arctan2(rot[1][0]/np.cos(B2), rot[0][0]/np.cos(B2))  
    return {"A1": deg*A1, "A2": deg*A2, "B1": deg*B1, "B2": deg*B2, "C1": deg*C1, "C2": deg*C2 }

def poly(x0, x1, x2, y0, y1, y2):
    A = np.array([[x0**2,x0,1],[x1**2,x1,1],[x2**2,x2,1]])
    B = np.array([y0,y1,y2])
    C = np.linalg.solve(A,B)   
    return C[0], C[1], C[2]

# only for second-degree polynomial function
@dataclass
class VectorFunction:
    x: sym.core.add.Add
    y: sym.core.add.Add
    z: sym.core.add.Add
    tangent: sym.core.add.Add = None
    norm: sym.core.add.Add = None
    binorm: sym.core.add.Add = None
    
    def __post_init__(self):
        t = sym.symbols('t')
        vec = [self.x, self.y, self.z]
        dvec = [sym.diff(coor_fun, t) for coor_fun in vec]
        mdvec = sym.sqrt(dvec[0]**2+dvec[1]**2+dvec[2]**2)    
        self.tangent = [dcoor_fun/mdvec for dcoor_fun in dvec]
    
    def tangent_val(self, val: float, coeffs: list) -> list:
        a0,a1,a2,t = sym.symbols('a0 a1 a2 t')
        return np.array([float(coor_fun.subs([(t,val),(a0,coeffs[0]),(a1,coeffs[1]),(a2,coeffs[2])])) for coor_fun in self.tangent])

