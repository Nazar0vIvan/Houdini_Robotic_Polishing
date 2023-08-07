# -*- coding: utf-8 -*-

import numpy as np
from math import * 

class TCP:
    def __init__(self,X,Y,Z,A,B,C,shoulder,elbow):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.A = A
        self.B = B
        self.C = C
        self.shoulder = shoulder
        self.shoulder = elbow

d = np.array([223.0, -101.5, 101.5, -660.0, 0.0, -80])
a = np.array ([150.0, -610.0, -20.0, 0.0, 0.0, 0.0])
alpha = np.array ([pi/2, 0, -pi/2, pi/2, -pi/2, pi])

def rotate(angle, axis):
    if (axis == "x"):
        rot = np.array([[1.0,0.0,0.0],[0.0,np.cos(angle),-np.sin(angle)],[0.0,sin(angle),np.cos(angle)]])        
    elif (axis == "y"):
        rot = np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-sin(angle),0,cos(angle)]])
    elif (axis == "z"):
        rot = np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])       
    rot[np.absolute(rot)<=0.0001] = 0.0
    return rot

def eulerAnglesToRatationMatrix(phi = 0, theta = 0, psi = 0):
    psi = psi*pi/180.0
    theta = theta*pi/180.0
    phi = phi*pi/180.0
    result = np.dot(np.dot(rotate(phi,"z"),rotate(theta,"y")),rotate(psi,"x"))
    result[np.absolute(result)<=0.0001]=0
    return result

# Denavit–Hartenberg parameters
d = np.array([223.0, -101.5, 101.5, -660.0, 0.0, -80])
a = np.array ([150.0, -610.0, -20.0, 0.0, 0.0, 0.0])
alpha = np.array ([pi/2, 0, -pi/2, pi/2, -pi/2, pi])   

# Direct Kinematics
def solveDK(nodeBase,nodeFlange,q2,q3):
    ABW_hou = nodeBase.worldTransform().transposed().asTupleOfTuples()
    ABW = np.asarray(ABW_hou)
    ABW[np.absolute(ABW)<=0.0001]=0
    
    A6W_hou = nodeFlange.worldTransform().transposed().asTupleOfTuples()
    A6W = np.asarray(A6W_hou)
    A6W[np.absolute(A6W)<=0.0001]=0
    
    A6B = np.dot(np.linalg.inv(ABW),A6W)    
    A6B[np.absolute(A6B)<=0.0001]=0

    #theta - pitch
    pitch = np.arctan2(-A6B[2,0],np.sqrt(A6B[2,1]**2+A6B[2,2]**2)) 
    
    #phi - roll
#    if ((pitch == pi/2) or (pitch == -pi/2)):
#        roll = 0
#    else:
    roll = np.arctan2(A6B[1,0]/np.cos(pitch),A6B[0,0]/np.cos(pitch))
        
    #psi - yaw
#    if (pitch == pi/2):
    yaw = np.arctan2(A6B[2,1]/np.cos(pitch),A6B[2,2]/np.cos(pitch))
        
    shoulder = "DOWN" if q2 > 0 else "UP"
    elbow =    "DOWN" if q3 >= 1.735 else "UP"           

    X = A6B[0,3]
    Y = A6B[1,3]
    Z = A6B[2,3]
    A = roll*180.0/pi
    B = pitch*180.0/pi 
    C = yaw*180.0/pi
 
    return {"TCP" : TCP(X,Y,Z,A,B,C,shoulder,elbow), "A6B" : A6B}

# Inverse Kinematics     
def solveIK(tcp, nodeBase, nodeRotatingColumn, nodeLinkArm, nodeWrist1):

    o6B = np.array([tcp.X,tcp.Y,tcp.Z,1])
    R6B = eulerAnglesToRatationMatrix(tcp.A,tcp.B,tcp.C)    
    A6B = np.column_stack((np.row_stack((R6B,np.zeros(3))),o6B))
    A6B[abs(A6B)<=0.0001] = 0
    
#    A6W_hou = nodeFlange.worldTransform().transposed().asTupleOfTuples()
#    A6W = np.asarray(A6W_hou)
#    A6W[np.absolute(A6W)<=0.0001]=0
#    A6B = np.dot(np.linalg.inv(ABW),A6W)
       
    # INVERSE POSITION #

    shoulder = 1.0 if tcp.shoulder == "DOWN" else -1.0
    elbow    = 1.0 if tcp.elbow == "DOWN" else -1.0
    
    # joint №1 #
    o6B = A6B[0:3,3]   
    o4B = o6B - abs(d[5])*A6B[0:3,2]
    q1 = np.arctan2(o4B[1],o4B[0])*180.0/pi
    
    # joint №2 #
    phi = np.arctan(abs(d[3])/abs(a[2]))*180.0/pi
    l = np.hypot(abs(d[3]),abs(a[2])) 
    
    A0B_hou = nodeRotatingColumn.localTransform().transposed().asTupleOfTuples()
    A0B = np.asarray(A0B_hou)
    A0B[np.absolute(A0B)<=0.0001]=0
    
    A60 = np.dot(np.linalg.inv(A0B),A6B)
    A60[np.absolute(A60)<=0.0001]=0
    
    A10_hou = nodeLinkArm.localTransform().transposed().asTupleOfTuples()
    A10 = np.asarray(A10_hou)
    A10[np.absolute(A10)<=0.0001]=0
            
    o60 = A60[0:3,3]    
    o40 = o60 - abs(d[5])*A60[0:3,2]    
    o10 = A10[0:3,3]
    o41 = o40 - o10
    o41m = np.hypot(o41[0],o41[2])
    
    gamma1 = np.arctan2(o41[2],o41[0])*180.0/pi
    gamma2 = np.arccos((abs(a[1])**2+o41m**2-l**2)/(2*abs(a[1])*o41m))*180.0/pi
    
    if ((shoulder == -1) and (elbow == 1) and (o41[2] < 0) and (o41[0] < 0)):
        q2 = -gamma1 - gamma2 - 360.0
    elif ((shoulder == -1) and (elbow == -1) and (o41[2] < 0) and (o41[0] < 0)):
        q2 = -gamma1 + gamma2 - 360.0
    else:
        q2 = -gamma1 - elbow*gamma2 
       
    # joint №3 #      
    beta = np.arccos((abs(a[1])**2+l**2-o41m**2)/(2*abs(a[1])*l))*180.0/pi    
    q3 = elbow*(180.0 - elbow*phi - beta + elbow*90.0)
            
    # INVERSE ORIENTATION #

    ABW_hou = nodeBase.worldTransform().transposed().asTupleOfTuples()
    ABW = np.asarray(ABW_hou)
    ABW[np.absolute(ABW)<=0.0001]=0

    R6B = A6B[0:3,0:3]
    A3W_hou = nodeWrist1.worldTransform().transposed().asTupleOfTuples()
    A3W = np.asarray(A3W_hou)
    A3W[np.absolute(A3W)<=0.0001]=0
    A3B = np.dot(np.linalg.inv(ABW),A3W)
    A36 = np.dot(np.linalg.inv(A3B),A6B)
    
#    A36 = np.dot(np.linalg.inv(A3W),A6W)
    
    o36 = A36[0:3,3]

#    R36 = A36[0:3,0:3]
#    R36[np.absolute(R36)<=0.001]=0
#    
#    A2W_hou = IN_LINE_WRIST.worldTransform().transposed().asTupleOfTuples()
#    A2W = np.asarray(A2W_hou)
#    A2W[np.absolute(A2W)<=0.0001]=0
#    
#    A32 = np.dot(np.linalg.inv(A2W),A3W)

    R1 = np.dot(np.dot(rotate(q1*pi/180.0,"z"),rotate(alpha[0],"y")),rotate(-alpha[0],"x"))
    R2 = np.dot(rotate(q2*pi/180.0,"z"),rotate(alpha[1],"x"))    
    R3 = np.dot(rotate(q3*pi/180.0,"z"),rotate(alpha[2],"x"))

    R3B = np.dot(np.dot(R1,R2),R3)    
    R36 = np.dot(np.transpose(R3B),R6B)
                   
#    # joint №4 #
#    q4 = -np.arctan2(A32[2,0],A32[0,0])*180.0/pi
    
    # joint №5 #       
    q5 = np.sign(o36[0])*np.arccos(-R36[2,2])*180.0/pi
    
    # joint №6 #
    if ((q5 >= 0.0) and (q5 < 180.0)):
        q6 = np.arctan2(R36[2,1],R36[2,0])*180.0/pi
        q4 = np.arctan2(R36[1,2],R36[0,2])*180.0/pi
    else:
        q6 = np.arctan2(-R36[2,1],-R36[2,0])*180.0/pi
        q4 = np.arctan2(-R36[1,2],-R36[0,2])*180.0/pi

    return {"q" : [0, q1,q2,q3,q4,q5,q6], "A6B" : A6B}

