import hou
import numpy as np
import sympy as sym
from math import *
import sys

from abc import ABCMeta, abstractmethod
from enum import Enum

deg = 180/pi

class TCP:
    def __init__(self,X,Y,Z,A,B,C,shoulder,elbow):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.A = A
        self.B = B
        self.C = C
        self.shoulder = shoulder
        self.elbow = elbow

d = np.array([223.0, -101.5, 101.5, -660.0, 0.0, -80])
a = np.array ([150.0, -610.0, -20.0, 0.0, 0.0, 0.0])
alpha = np.array ([pi/2, 0, -pi/2, pi/2, -pi/2, pi]) 

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
    
# Denavit–Hartenberg parameters
d =     [0.,        450.,    0.,   0., -660., 0., 80.]
a =     [0.,          0.,  150., 610.,   20., 0.,  0.]
alpha = [0.,          0., -pi/2.,   0.0, pi/2., -pi/2., -pi/2]
q0 =    [0.,          0.,     0., -pi/2,    0.,     0.,  pi]

def dhMatrix(q,i):
    qf = q + q0[i];
    mat = np.array([[cos(qf),               -sin(qf),                0,              a[i]],
                    [cos(alpha[i])*sin(qf),  cos(alpha[i])*cos(qf), -sin(alpha[i]), -d[i]*sin(alpha[i])],
                    [sin(alpha[i])*sin(qf),  sin(alpha[i])*cos(qf),  cos(alpha[i]),  d[i]*cos(alpha[i])],
                    [0,                      0,                      0,              1]]);
    return mat
    
mm = np.array([[1,2,3,4],[2,2,2,2],[3,3,3,3],[4,4,4,4]])

AF0 = np.array([[ 0.0272, 0.9957,  0.0888, 698.1308],
                [ 0.0571, 0.0871, -0.9946, 198.7530],
                [-0.9980, 0.0321, -0.0545, 643.1010],
                [ 0,      0,       0,        1.0000]])

shoulder = -1
elbow = 1
                
# Inverse Kinematics
def solveIK(AF0, shoulder, elbow):
# q1
    o600 = AF0[0:3,[3]]
    o400 = o600 - abs(d[6])*AF0[0:3,[2]]
    q1 = np.arctan2(o400[1][0],o400[0][0])
# q2
    phi = np.arctan(abs(d[4])/abs(a[4]))
    L = hypot(abs(d[4]),abs(a[4]))
    
    A10 = dhMatrix(-q1,1)
    o400r = np.dot(A10[0:3,0:3],o400) # o400r - вектор, проведенный в начало 4-ой СК из начала 0-ой СК и заданный в 0-ой СК, но повернутой на угол q1
    
    o421s = o400r - np.array([[150.],[0.],[450.]]) # o421s - вектор, проведенный в начало 4-ой СК из начала 2-ой СК и заданный в 1-ой СК, но сдвинутой в начало 2-ой СК
    r42m = hypot(o421s[0],o421s[2]) # o421s - модуль вектора o421
    
    gamma1 = np.arctan2(o421s[2][0],o421s[0][0])
    gamma2 = np.arccos((abs(a[3])**2 + r42m**2 - L**2)/(2*abs(a[3])*r42m))
    
    shoulder = 1.0 if shoulder == "DOWN" else -1.0
    elbow    = 1.0 if elbow == "DOWN" else -1.0
    
    if ((shoulder == -1) and (elbow == 1) and (o421s[2] < 0) and (o421s[0] < 0)):
        q2 = -gamma1 - gamma2 - 2*pi
    elif ((shoulder == -1) and (elbow == -1) and (o421s[2] < 0) and (o421s[0] < 0)):
        q2 = -gamma1 + gamma2 - 2*pi
    else:
        q2 = -gamma1 - elbow*gamma2

# q3
    beta = np.arccos((abs(a[3])**2+L**2-r42m**2)/(2*abs(a[3])*L))
    q3 = elbow*(pi - elbow*phi - beta + elbow*pi/2)

# q4, q5, q6
    A30 = np.dot(dhMatrix(q1,1),np.dot(dhMatrix(q2,2),dhMatrix(q3,3)))
    A63 = np.dot(np.linalg.inv(A30),AF0)

    q4 = np.array([0.,0.]); q5 = np.array([0.,0.]); q6 = np.array([0.,0.])
    
    # q5
    q5[0] = np.arccos(A63[1,2])
    q5[1] = -q5[0]
    
    # q4
    q4[0] = np.arctan2(-A63[2,2]/sin(q5[0]),-A63[0,2]/sin(q5[0]))
    q4[1] = np.arctan2(-A63[2,2]/sin(q5[1]),-A63[0,2]/sin(q5[1]))

    # q6
    q6[0] = np.arctan2(A63[1,1]/sin(q5[0]),-A63[1,0]/sin(q5[0]))
    q6[1] = np.arctan2(A63[1,1]/sin(q5[1]),-A63[1,0]/sin(q5[1]))
    
    shoulder = "DOWN" if q2 > 0 else "UP"
    elbow =    "DOWN" if q3 >= 1.735 else "UP"  

    return{ 'qt': [q1*deg,q2*deg,q3*deg], 
            'qr': [q4*deg,q5*deg,-q6*deg],
            'shoulder': shoulder,
            'elbow': elbow }
    
'''
def solveIK(tcp, nodeBase, nodeRotatingColumn, nodeLinkArm, nodeWrist1):

    o6B = np.array([tcp.X,tcp.Y,tcp.Z,1])
    R6B = eulerAnglesToRotationMatrix(tcp.A,tcp.B,tcp.C)    
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

    R1 = np.dot(np.dot(rotationMatrix(q1*pi/180.0,"z"),rotationMatrix(alpha[0],"y")),rotationMatrix(-alpha[0],"x"))
    R2 = np.dot(rotationMatrix(q2*pi/180.0,"z"),rotationMatrix(alpha[1],"x"))    
    R3 = np.dot(rotationMatrix(q3*pi/180.0,"z"),rotationMatrix(alpha[2],"x"))

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
    
    shoulder = "DOWN" if q2 > 0 else "UP"
    elbow =    "DOWN" if q3 >= 1.735 else "UP"  

    return {"q" : [0, q1,q2,q3,q4,q5,q6], "A6B" : A6B, "shoulder" : shoulder, "elbow" : elbow}
'''

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
# ----------------------------------------------------------------------------------
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

### PLANE ###
def points2plane(x,y,z):
    n = len(x)

    U = np.zeros((3,3))
    V = np.zeros(3)
    
    U[0,0] = np.sum(x**2)
    U[0,1] = np.sum(x*y)
    U[0,2] = np.sum(x)
    
    U[1,0] = U[0,1]
    U[1,1] = np.sum(y**2)
    U[1,2] = np.sum(y)
    
    U[2,0] = U[0,2]
    U[2,1] = U[1,2]
    U[2,2] = n
    
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
    tau = []; bnorm = []
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

def drawPrimitivePoints(geo, primitivePoints):
    for i in range(len(primitivePoints.points)):
        point = geo.createPoint()
        point.setPosition((primitivePoints.points[i].x, primitivePoints.points[i].y, primitivePoints.points[i].z))
        
def drawRectangle(geo, p1, p2, p3, p4): 
    poly = geo.createPolygon()
    poly.addVertex(p1)
    poly.addVertex(p2)
    poly.addVertex(p3)
    poly.addVertex(p4)

### FRENE FRAME: END ###

def drawPoint(geo, point):
    p = geo.createPoint()
    p.setPosition((point.x, point.y, point.z))
   
np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)