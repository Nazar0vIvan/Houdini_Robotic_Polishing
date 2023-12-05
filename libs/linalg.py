import numpy as np
from math import pi, sin, cos

def rotationMatrix3x3(angle, axis):
    if (axis == "x"):
        rot = np.array([[1.,0.,0.],[0.,np.cos(angle),-np.sin(angle)],[0.,sin(angle),np.cos(angle)]])        
    elif (axis == "y"):
        rot = np.array([[np.cos(angle),0.,np.sin(angle)],[0.,1.,0.],[-sin(angle),0.,cos(angle)]])
    elif (axis == "z"):
        rot = np.array([[np.cos(angle),-np.sin(angle),0.],[np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])       
    rot[np.absolute(rot)<=0.0001] = 0.
    return rot
    
def rotationMatrix4x4(angle, axis):
    if (axis == "x"):
        rot = np.array([[1.,0.,0.,0.],[0.,np.cos(angle),-np.sin(angle), 0.],[0.,sin(angle),np.cos(angle),0.],[0.,0.,0.,1.]])        
    elif (axis == "y"):
        rot = np.array([[np.cos(angle),0.,np.sin(angle),0.],[0.,1.,0.,0.],[-sin(angle),0.,cos(angle),0.],[0.,0.,0.,1.]])
    elif (axis == "z"):
        rot = np.array([[np.cos(angle),-np.sin(angle),0.,0.],[np.sin(angle),np.cos(angle),0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])       
    rot[np.absolute(rot)<=0.0001] = 0.0
    return rot

def translationMatrix(dx,dy,dz):
    return np.array([[1.,0.,0.,dx],[0.,1.,0.,dy],[0.,0.,1.,dz],[0.,0.,0.,1.]])

def euler2rot(A = 0, B = 0, C = 0):
    A = A*pi/180.0; B = B*pi/180.0; C = C*pi/180.0
    result = np.dot(np.dot(rotationMatrix3x3(A, "z"), rotationMatrix3x3(B, "y")), rotationMatrix3x3(C, "x"))
    result[np.absolute(result)<=0.0001]=0
    return result

def rot2euler(rot, is_deg = False):
    deg = 180/pi if is_deg else 1
    B1 = -np.arcsin(rot[2][0])
    B2 = pi + np.arcsin(rot[2][0])
    C1 = np.arctan2(rot[2][1]/np.cos(B1), rot[2][2]/np.cos(B1))
    C2 = np.arctan2(rot[2][1]/np.cos(B2), rot[2][2]/np.cos(B2))
    A1 = np.arctan2(rot[1][0]/np.cos(B1), rot[0][0]/np.cos(B1))
    A2 = np.arctan2(rot[1][0]/np.cos(B2), rot[0][0]/np.cos(B2))  
    return {"A1": deg*A1, "A2": deg*A2, "B1": deg*B1, "B2": deg*B2, "C1": deg*C1, "C2": deg*C2 }

