from dataclasses import dataclass
import numpy as np

@dataclass
class Frene:
    t: np.ndarray
    b: np.ndarray
    n: np.ndarray
    p: np.ndarray
    transf: np.ndarray = None

    def __post_init__(self):
        self.transf = np.column_stack((np.append(self.t, 0.), np.append(self.b, 0.), np.append(self.n, 0.), np.append(self.p, 1.)))        

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
