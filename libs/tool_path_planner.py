import numpy as np
from math import pi, sin, cos
from dataclasses import dataclass
from math import *
from linalg import *

# CLASSES

@dataclass
class Roll:
    r: float
    h: float
    
class Airfoil:
    def __init__(self, x: list, y: list, z: list, npf: int, npc: int, npe: int, rz: float):
        self.npf = npf
        self.npc = npc
        self.npe = npe
        self.rz = rz
        npt = 2*npc + 2*npe
        self.npt = npt
        self.geo = {}
        for i in range(self.npf):
            self.geo[i] = []
            for j in range(npt):
                self.geo[i].append(np.array([x[i*npt+j], y[i*npt+j], z[i*npt+j]]))              
        self.cx_length = sum(np.linalg.norm(self.geo[0][i+1] - self.geo[0][i]) for i in range(self.npc-1)) 
        self.airfoil_height = np.linalg.norm(self.geo[0][0] - self.geo[1][0])

@dataclass
class Frene:
    t: np.ndarray
    b: np.ndarray
    n: np.ndarray
    p: np.ndarray
    transform: np.ndarray = None
    def __post_init__(self):
        self.transform = np.column_stack((np.append(self.t,0.), np.append(self.b,0.), np.append(self.n,0.), np.append(self.p,1.)))        
        
class RuledSurf:
    def __init__(self, a1: np.ndarray, a2: np.ndarray, b1: np.ndarray, b2: np.ndarray):
        self.a1 = a1
        self.a2 = a2
        self.b1 = b1
        self.b2 = b2
    def val(self, u: float, v: float) -> np.array:
        return (1-u)*(1-v)*self.a1 + u*(1-v)*self.a2 + (1-u)*v*self.b1 + u*v*self.b2      
    def du(self, v: float) -> np.array:
        du = (self.a1 - self.a2)*(v-1) - (self.b1 - self.b2)*v
        return du/np.linalg.norm(du)            
    def dv(self, u: float) -> np.array:
        dv = (self.a1 - self.b1)*(u-1) - (self.a2 - self.b2)*u
        return dv/np.linalg.norm(dv)
    def get_frene(self, u: float, v: float) -> Frene:
        pt = self.val(u,v)
        du = self.du(v)
        b = -self.dv(u)
        n = np.cross(du, b)
        t = np.cross(b, n)
        return Frene(t, b, n, pt)
    def get_ufrenes(self, v: float, nu: int):
        ufrenes = []
        step_u = 1./nu
        u = 0.
        for iu in range(nu):
            ufrenes.append(self.get_frene(u,v))
            u = u + step_u
        return ufrenes
        
class ToolPathPlanner: 
    def __init__(self, fps, path = None):
        self.__path = [] if path == None else path
        self.__fps = fps
        
    @property
    def path(self):
        return self.__path
        
    def add_lead_in(self, cut_path: list, vel: float, lead_in_dist: float, angle = 0.0) -> None:
        p1 = np.array(cut_path[0])
        p2 = np.array(cut_path[1])     
        v21 = p1[0:3] - p2[0:3]
        v13 = lead_in_dist/np.linalg.norm(v21)*v21
        p3 = np.concatenate((p1[0:3] + v13, p1[3:6]))
        
        length = np.linalg.norm(p1[0:3] - p3[0:3])
        nframes = floor(length/vel*self.__fps)
        
        self.__path.extend(np.linspace(p3.tolist(), p1.tolist(), nframes+2, endpoint=False))    
        
    def add_lead_out(self, cut_path: list, vel: float, lead_out_dist: float, angle = 0.0) -> None:
        p1 = np.array(cut_path[-1])
        p2 = np.array(cut_path[-2])
        v21 = p1[0:3] - p2[0:3]
        v13 = lead_out_dist/np.linalg.norm(v21)*v21
        p3 = np.concatenate((p1[0:3] + v13, p1[3:6]))
        
        length = np.linalg.norm(p1[0:3] - p3[0:3])
        nframes = floor(length/vel*self.__fps)
        
        self.__path.extend(np.linspace(p1.tolist(), p3.tolist(), nframes+2, endpoint=False).tolist())
    
    def add_lin(self, p1: list, p2: list, vel: float, is_const_orien = False) -> None:
        length = np.linalg.norm(p2[0:3] - p1[0:3])
        nframes = length/vel*self.__fps
        if(is_const_orien):
            lin_path = np.linspace(p1, p2, nframes+2, endpoint=False).tolist()
        else:
            lin_path = np.linspace(p1[0:3], p2[0:3], nframes+2, endpoint=False)
            lin_path = [p.extend(p1[3:6]) for p in lin_path]
        self.__path.extend(lin_path)
    def add_path(self, path_to_add: list) -> None:
        self.__path.extend(path_to_add)

# FUNCTIONS

def get_stripe_tool_path(v: float, nu: int, airfoil: Airfoil, roll: Roll) -> list:    
    Tr_C_i = get_rhs_transf(roll)
    frenes = []
    for i in range(airfoil.npc-1): #airfoil.npc-1
        ruled_surf = get_ruled_surf(airfoil, i)
        frenes.extend(ruled_surf.get_ufrenes(v, nu))
    return frenes2path(frenes, Tr_C_i) # Tr_B_0

def get_rhs_transf(roll: Roll, is_cls = False) -> np.ndarray:
    if (is_cls):
        Tr_S_i = np.array([[1.,0.,0.,0.], [0.,0.,1.,0.], [0.,-1.,0.,0.], [0.,0.,0.,1.]])
        Tr_C_S = np.array([[1.,0.,0.,0.], [0.,1.,0.,-roll.r], [ 0.,0.,1.,0.], [0.,0.,0.,1.]])
    else:
        Tr_S_i = np.array([[0.,0.,1.,0.], [0.,1.,0.,0.], [-1.,0.,0.,0.], [0.,0.,0.,1.]])
        Tr_C_S = np.array([[1.,0.,0.,-roll.r], [0.,1.,0.,roll.h/2], [ 0.,0.,1.,0.], [0.,0.,0.,1.]])
    return np.dot(Tr_S_i, Tr_C_S)

def get_ruled_surf(airfoil: dict, j: int, i = 0) -> RuledSurf:
    a1 = airfoil.geo[i][j];   a2 = airfoil.geo[i][j+1]
    b1 = airfoil.geo[i+1][j]; b2 = airfoil.geo[i+1][j+1]
    return RuledSurf(a1,a2,b1,b2)

def frenes2path(frenes: list,  rhs_transf: np.ndarray) -> dict:
    tool_path = []
    transforms = []
    for frene in frenes:
        #Tr_C_0 = np.dot(lhs_transf, np.dot(frene.transform, rhs_transf))
        Tr_C_0 = np.dot(frene.transform, rhs_transf)
        transforms.append(Tr_C_0)
        euler = rot2euler(Tr_C_0[0:3,0:3], True)
        tool_path.append([*Tr_C_0[0:3,3], euler["C1"], euler["B1"], euler["A1"]])
    return tool_path, transforms