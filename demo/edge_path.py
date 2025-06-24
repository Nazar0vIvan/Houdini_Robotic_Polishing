import numpy as np
import json
from linalg import *
from diffgeom import *
from math import floor

def get_frene(ptc: np.array, pt: np.array):
    n = pt - ptc
    b = np.array([0.,0.,1.])
    n = n/np.linalg.norm(n) # normal
    t = np.cross(b,n)
    return Frene(np.round(t,3), np.round(b,3), np.round(n,3), pt)

def get_frenes(blade: dict, edge_type: str) -> list:
    frenes = []
    npf = len(blade)
    npe = len(blade[0]["re"])
    for i in range(npf):
        ptc = blade[i][f"{edge_type}_edge"][0:3]
        pt = blade[i][edge_type][floor(npe/2)]
        frenes.append(get_frene(ptc, pt))
    return frenes

################# BEGIN #################

with open("99.01.25.242.json", 'r') as file:
    blade = json.load(file)

lead_dist = 30
pe_count = len(blade[0]["re"])

AiL = np.array([[0.,1.,0.,0.],[1.,0.,0.,0.],[0.,0.,-1.,0.],[0.,0.,0.,1.]]) # WHEEL

re_frenes = get_frenes(blade, "re")
ABTs = []
path = []
#i = 0; j = 0
for pf_frene in re_frenes:
    pf_path = []
    #j = 0
    for frene in pf_frene:
        ABT = np.dot(AiL, np.linalg.inv(frene.transf))
        ABTs.append(ABT)
        euler = rot2euler(ABT, True)
        A = euler.get('A1')
        B = euler.get('B1')
        C = euler.get('C1')
        X = ABT[0][3]
        Y = ABT[1][3]
        Z = ABT[2][3]
        pf_path.append([X,Y,Z,A,B,C])
        #print(f"{i}|{j}: ", np.around(trajectory[-1], decimals=2))
        #j = j + 1
    #i = i + 1
    path.append(pf_path)

rounded_path = np.round(np.array(path), 3)
