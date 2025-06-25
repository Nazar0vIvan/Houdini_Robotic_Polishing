import numpy as np
import json
from linalg import *
from diffgeom import *
from math import floor

def get_lead_point(ptc: np.array, pt: np.array, dist: float) -> np.ndarray:
    #ptc - edge circ center
    #pt - point of the circ edge
    v21 = pt - ptc
    v13 = dist/np.linalg.norm(v21)*v21
    return np.round(pt + v13, 3)

def get_frene(ptc: np.array, pt: np.array):
    n = pt - ptc
    b = np.array([0.,0.,1.])
    n = n/np.linalg.norm(n) # normal
    t = np.cross(b,n)
    return Frene(np.round(t,3), np.round(b,3), np.round(n,3), pt)

def get_frenes(blade: dict, edge_type: str, pf_num_begin: int, pf_num_end: int, lead_dist: float) -> list:
    frenes = [Frene()]
    npe = len(blade[0]["re"])
    for i in range(pf_num_begin, pf_num_end+1):
        ptc = blade[i][f"{edge_type}_edge"][0:3]
        pt = blade[i][edge_type][floor(npe/2)]
        frenes.append(get_frene(np.array(ptc), np.array(pt)))
    lead_in = get_lead_point(np.array(blade[pf_num_begin][f"{edge_type}_edge"][0:3]), blade[pf_num_begin][edge_type][floor(npe/2)], lead_dist)
    lead_out = get_lead_point(np.array(blade[pf_num_end][f"{edge_type}_edge"][0:3]), blade[pf_num_end][edge_type][floor(npe/2)], lead_dist)
    frenes[0] = Frene(frenes[1].t, frenes[1].b, frenes[1].n, lead_in)
    frenes.append(Frene(frenes[-1].t, frenes[-1].b, frenes[-1].n, lead_out))
    return frenes

def generate_edge_dat(program_name: str, path_name_prefix: str, path: list) -> None:
    n = len(path)
    with open(f"{program_name}.dat", 'w') as file:
        file.write(f"DECL POS {path_name_prefix}[{n}]\n\n")
        for i, pos in enumerate(path):
            file.write(f"{path_name_prefix}[{i+1}] = {{X {pos[0]}, Y {pos[1]}, Z {pos[2]}, A {pos[3]}, B {pos[4]}, C {pos[5]}}}\n")

def generate_edge_src(program_name: str, path_name_prefix: str, path: list) -> None:
    n = len(path)
    with open(f"{program_name}.src", 'w') as file:
        file.write("$BASE = BASE_DATA[6]\n")
        file.write("$TOOL = TOOL_DATA[2]\n")
        file.write("LIN_REL {X -200, Y -500}\n")
        file.write(f"PTP {path_name_prefix}[{1}]\n")
        file.write("SPLINE WITH $VEL= {CP 0.007, ORI1 0.005, ORI2 0.005}\n")
        file.write(f"  SLIN {path_name_prefix}[{2}]\n")
        for i in range(3,n):
            file.write(f"  SPL {path_name_prefix}[{i}]\n")
        file.write(f"  SLIN {path_name_prefix}[{n}] WITH $VEL.CP = 0.01\n")
        file.write("ENDSPLINE\n\n")
        file.write("$VEL.CP = 0.07\n")
        file.write("LIN_REL {Z 400}")

################# BEGIN #################

with open("99.01.25.242.json", 'r') as file:
    blade = json.load(file)

lead_dist = 20
pf_num_begin = 3
pf_num_end = 18
AiL = np.array([[0.,1.,0.,0.],[1.,0.,0.,0.],[0.,0.,-1.,0.],[0.,0.,0.,1.]]) # WHEEL

re_frenes = get_frenes(blade, "re", pf_num_begin, pf_num_end, lead_dist)
ABTs = []
re_path = []
#i = 0; j = 0
for frene in re_frenes:
    #j = 0
    ABT = np.dot(AiL, np.linalg.inv(frene.transf))
    ABTs.append(ABT)
    euler = rot2euler(ABT, True)
    A = euler.get('A1')
    B = euler.get('B1')
    C = euler.get('C1')
    X = ABT[0][3]
    Y = ABT[1][3]
    Z = ABT[2][3]
    re_path.append([X,Y,Z,A,B,C])
    #print(f"{i}|{j}: ", np.around(trajectory[-1], decimals=2))
    #j = j + 1
    #i = i + 1

re_rounded_path = np.round(np.array(re_path), 3)

generate_edge_dat("re_wheel", "RE", re_rounded_path)
generate_edge_src("re_wheel", "RE", re_rounded_path)


