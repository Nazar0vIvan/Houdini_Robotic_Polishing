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

def get_frene(pt0: np.ndarray, u1: np.ndarray, u2: np.ndarray, v1: np.ndarray):
    circ = points2circ(pt0, u1, u2)
    ptc = circ.center
    n = pt0 - ptc
    n = n/np.linalg.norm(n)
    b = pt0 - v1
    b = b/np.linalg.norm(b) # normal
    t = np.cross(b,n)
    return Frene(np.round(t,3), np.round(b,3), np.round(n,3), pt0)

def get_frenes(blade: dict, edge_type: str, pf_num_begin: int, pf_num_end: int, lead_dist: float) -> list:
    npe = len(blade[0]["re"])
    frenes = []
    for i in range(pf_num_begin, pf_num_end+1):
        u1 = blade[i][edge_type][0]
        pt0 = blade[i][edge_type][floor(npe/2)]
        u2 = blade[i][edge_type][-1]
        v1 = blade[i+1][edge_type][floor(npe/2)]
        frenes.append(get_frene(np.array(pt0), np.array(u1), np.array(u2), np.array(v1)))
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
pf_num_begin = 0
pf_num_end = 70
AiL = np.array([[0.,1.,0.,0.],[1.,0.,0.,0.],[0.,0.,-1.,0.],[0.,0.,0.,1.]]) # WHEEL

re_frenes = get_frenes(blade, "re", pf_num_begin, pf_num_end, lead_dist)
ABTs = []
re_path = []
re_path = [[0.,0.,0.,0.,0.,0.]]
for frene in re_frenes:
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

first = re_path[1]
re_path[0] = [first[0], first[1], first[2] + 30, first[3], first[4], first[5]]
re_rounded_path = np.round(np.array(re_path), 3)

generate_edge_dat("re_wheel", "re", re_rounded_path)
generate_edge_src("re_wheel", "re", re_rounded_path)


