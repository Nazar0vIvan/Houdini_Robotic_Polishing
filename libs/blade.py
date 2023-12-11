from dataclasses import dataclass
import numpy as np
from diffgeom import Frene

# npf - number of profiles
# npc - number of points on the convex/concave
# npe - number of points on the edge
# npr - number of points on radius sector

@dataclass  
class Blade:
    convex: dict 
    radcx: dict
    npf: int
    npc: int
    npe: int
    npr: int
    airheight: float = 0.0

    def __post_init__(self):
        self.airheight = sum(np.linalg.norm(self.convex[pfnum][0] - self.convex[pfnum+1][0]) for pfnum in range(self.npf-1))
    
    def pf_cx_length(self, pfnum: int):
        return sum(np.linalg.norm(self.convex[pfnum][ptnum+1] - self.convex[pfnum][ptnum]) for ptnum in range(self.npc-1))
    
    def pf_radcx_length(self, pfnum: int):
        return sum(np.linalg.norm(self.radcx[pfnum][ptnum+1] - self.radcx[pfnum][ptnum]) for ptnum in range(self.npc-1)) 