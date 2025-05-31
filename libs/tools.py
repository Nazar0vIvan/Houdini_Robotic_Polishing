import numpy as np
from dataclasses import dataclass
from math import tan, radians
from linalg import *


# {T}  - tool frame
# {TS} - frame, attached to the point on the tool cutting surface
# {BS} - frame, attached to the point on the blank surface

class Tool:
    def __init__(self):
        self._transf_surf2blank = np.eye(4, dtype=float) # {TS} -> {BS}
        self.transf_orig2surf = np.eye(4, dtype=float) # {T} -> {TS}
        self.transf_orig2blank = np.eye(4, dtype=float) # {T} -> {BS}
        self.transf_orig2base = np.eye(4, dtype=float) # {T} -> {0}

    @property
    def transf_surf2blank(self) -> np.ndarray:
        return self._transf_surf2blank
        
    @transf_surf2blank.setter
    def transf_surf2blank(self, v: np.ndarray) -> None:
        self._transf_surf2blank = v
        self.transf_orig2blank = np.dot(self._transf_surf2blank, self.transf_orig2surf)

@dataclass
class W1FF1(Tool):
    D: float 
    B: float
    
    def __post_init__(self):
        super().__init__()
        self.transf_orig2surf = np.array([[1.,0.,0.,-self.D/2], [0.,1.,0.,0.], [0.,0.,1.,-self.B/2.], [0.,0.,0.,1.]])

@dataclass
class W5CHIX(Tool):
    D: float
    R: float
    H: float
    d: float
    D1: float
    D3: float
    DD: float = 0.0

    
    def __post_init__(self):
        super().__init__()
        alpha = 55.0
        beta = 25.0
        

        # full diameter DD calculation
        '''
        dD = (self.D3 - self.D1)/2
        dH = dD/tan(radians(90 - self.alpha))
        HH = dH + self.H
        left = HH*tan(radians(90 - self.alpha))/(tan(radians(90 - self.alpha)) + tan(radians(90 - self.beta)))
        dd = left*tan(radians(90 - self.beta))
        self.DD = self.D1 + 2*dd
        '''

@dataclass
class Roll(Tool):
    R: float
    H: float
    
    def __post_init__(self):
        super().__init__()
    
    '''
    def transf(self, is_cls = False) -> np.ndarray:
        if (is_cls):
            Tr_S_i = np.array([[1.,0.,0.,0.], [0.,0.,1.,0.], [0.,-1.,0.,0.], [0.,0.,0.,1.]])
            Tr_C_S = np.array([[1.,0.,0.,0.], [0.,1.,0.,-self.r], [ 0.,0.,1.,0.], [0.,0.,0.,1.]])
        else:
            Tr_S_i = np.array([[0.,0.,1.,0.], [0.,1.,0.,0.], [-1.,0.,0.,0.], [0.,0.,0.,1.]])
            Tr_C_S = np.array([[1.,0.,0.,-self.r], [0.,1.,0.,self.h/2], [ 0.,0.,1.,0.], [0.,0.,0.,1.]])
        return np.dot(Tr_S_i, Tr_C_S)
    '''


