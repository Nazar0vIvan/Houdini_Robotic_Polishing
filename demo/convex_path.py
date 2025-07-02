import numpy as np
import json
from linalg import *
from diffgeom import *

################# BELT {L} #################

# L - BELT: belt frame - t, b, n
def get_belt_frame(o: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray):
    [A, B, C, D, AA, BB, DD] = points2plane(x, y, z)

    n = np.array([A, B, C]) # belt z - n
    n = n / np.linalg.norm(n)

    x0 = np.array([1, 0, 0]) # X0 axis unit vector
    t = x0 - np.dot(x0, n) * n # belt x - t
    t = -t / np.linalg.norm(t)

    b = np.cross(n, t) # belt y - b

    transform = np.eye(4, dtype=np.float64)
    transform[:3, 0] = t
    transform[:3, 1] = b
    transform[:3, 2] = n
    transform[:3, 3] = o
    
    eulerA = np.degrees(np.arctan2(transform[1,0], transform[0,0]))
    eulerB = np.degrees(np.arcsin(-transform[2,0]))
    eulerC = np.degrees(np.arctan2(transform[2,1], transform[2,2]))
    frame = np.append(o, [eulerA, eulerB, eulerC])

    return [frame, transform]

################# BLADE ################# 

def get_lead_point(pt1: np.ndarray, pt2: np.ndarray, dist: float) -> np.ndarray:
    #pt1 - the last point of cx/cv
    #pt2 - the point before the last of cx/cv
    v = pt1 - pt2
    magnitude = np.linalg.norm(v)
    unit_v = v / magnitude
    return np.round(pt1 + unit_v * dist, 3)

def add_lead_points(blade: dict, lead_dist: float) -> dict:
    for pf in blade:
        pf["cx"][0] = get_lead_point(np.array(pf["cx"][1]), np.array(pf["cx"][2]), lead_dist)
        pf["cx"][-1] = get_lead_point(np.array(pf["cx"][-2]), np.array(pf["cx"][-3]), lead_dist)
    return blade


def get_frene(p0: np.array, u1: np.array, u2: np.array, v1: np.array) -> Frene:
    a0,a1,a2,t = sym.symbols('a0 a1 a2 t')
    xu = t; yu = a0*t**2+a1*t+a2; zu = u1[2];             
    vecfun_u = VectorFunction(xu, yu, zu)
    coeffs_u = poly(u1[0], p0[0], u2[0], u1[1], p0[1], u2[1]) 
    tanu = vecfun_u.tangent_val(p0[0], coeffs_u)
    if(tanu[0] < 0): tanu = -1*tanu
    tanv = v1 - p0
    tanv = tanv/np.linalg.norm(tanv)
    norm = np.cross(tanu, tanv)
    norm = norm/np.linalg.norm(norm)
    binorm = np.cross(norm, tanu)
    return Frene(tanu, binorm, norm, p0)

def get_cx_frenes(blade: dict, pf_num_begin: int, pf_num_end: int) -> list:
    frenes = []
    for i in range(pf_num_begin, pf_num_end):
        frenes.append(get_cx_airstripe_frenes(blade, i))
    return frenes

def get_cx_airstripe_frenes(blade: dict, pf_num: int) -> list:
    pf_frenes = [Frene()]
    npt = len(blade[pf_num]["cx"])
    for ptnum in range(1, npt-1):
        p0 = np.array(blade[pf_num]["cx"][ptnum])
        u1 = np.array(blade[pf_num]["cx"][ptnum-1])
        u2 = np.array(blade[pf_num]["cx"][ptnum+1])   
        v1 = np.array(blade[pf_num+1]["cx"][ptnum])
        pf_frenes.append(get_frene(p0, u1, u2, v1))
    pf_frenes[0] = Frene(pf_frenes[1].t, pf_frenes[1].b, pf_frenes[1].n, np.array(blade[pf_num]["cx"][0]))
    pf_frenes.append(Frene(pf_frenes[-1].t, pf_frenes[-1].b, pf_frenes[-1].n, np.array(blade[pf_num]["cx"][-1])))
    return pf_frenes

def generate_cx_dat(program_name: str, pf_prefix: str, path: list) -> None:
    n = len(path)
    m = len(path[0])
    with open(f"{program_name}.dat", 'w') as file:
        for i in range(n):
            file.write(f"\n\n; ### {pf_prefix}{i+1} ###\n\n")
            file.write(f"DECL POS {pf_prefix}{i+1}[{m}]\n\n")
            for j, pos in enumerate(path[i]):
                file.write(f"{pf_prefix}{i+1}[{j+1}] = {{X {pos[0]}, Y {pos[1]}, Z {pos[2]}, A {pos[3]}, B {pos[4]}, C {pos[5]}}}\n")
        
def generate_cx_src(program_name: str, pf_prefix: str, path: list) -> None:
    n = len(path)
    m = len(path[0])
    with open(f"{program_name}.src", 'w') as file:
        file.write("$BASE = BASE_DATA[3]\n")
        file.write("$TOOL = TOOL_DATA[2]\n")
        file.write("PTP {A4 90, A5 -90}")
        for i in range(n):
            file.write(f"\n\n; ### {pf_prefix}{i+1} ###\n")
            file.write(f"PTP {pf_prefix}{i+1}[{1}]\n")
            file.write("SPLINE WITH $VEL= {CP 0.007, ORI1 0.005, ORI2 0.005}\n")
            file.write(f"  SLIN {pf_prefix}{i+1}[{2}] WITH $VEL.CP = 0.01\n")
            for j in range(m-3):
                file.write(f"  SPL {pf_prefix}{i+1}[{j+3}]\n")
            file.write(f"  SPL {pf_prefix}{i+1}[{m}] WITH $VEL.CP = 0.01\n")
            file.write("ENDSPLINE\n\n")
            file.write("$VEL.CP = 0.07\n")
            file.write("LIN_REL {Z 400}") if (i == n-1) else file.write("LIN_REL {Z 20}")


def get_cx_path(cx_frenes: np.ndarray, AiL: np.ndarray):
    ABTs = []; path = []
    for pf_frene in cx_frenes:
        pf_path = []
        for frene in pf_frene:
            ABT = np.dot(AiL, np.linalg.inv(frene.transf))
            ABTs.append(ABT)     
            X, Y, Z = ABT[:3, 3]
            euler = rot2euler(ABT, True)
            A = euler['A1']; B = euler['B1']; C = euler['C1']
            pf_path.append([X, Y, Z, A, B, C])
        path.append(pf_path)
    return ABTs, path


#########################################
################# BEGIN #################
#########################################

# BELT_FRAME: L -> 0
oT = np.array([1009.15, -16.49, 623.81])
xT = np.array([996.14, 1010.89, 1010.89, 1023.99, 1014.15, 1014.15, 1004.89, 1004.89, 1009.15])
yT = np.array([-16.14, -29.24, 0.92, -16.14, -10.54, -22.95, -22.21, -10.51, -16.49])
zT = np.array([625.57, 623.52, 623.48, 622.35, 623.61, 622.86, 624.73, 624.40, 623.81])
[BELT_FRAME, AL0] = get_belt_frame(oT, xT, yT, zT)

np.set_printoptions(formatter={'float': '{:.3f}'.format})
print("BELT_FRAME:\n", "X: ", BELT_FRAME[0], "| Y: ", BELT_FRAME[1], "| Z: ", BELT_FRAME[2], "| A: ", BELT_FRAME[3], "| B: ", BELT_FRAME[4], "| C: ", BELT_FRAME[5])

with open("99.01.25.242.json", 'r') as file:
    blade = json.load(file)
# [
#   [Frene_11, Frene_12 ...]
#   [Frene_21, Frene_22 ...]
#   ...
#   [Frene_n1, Frene_n2 ...]
# ]

# B -> F
ABF = np.dot(translationMatrix([0.011, 0.047, 153.319]), rotationMatrix4x4(np.radians(-49), "z"))

# i -> T
#AiL = np.array([[-1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,-1.,0.],[0.,0.,0.,1.]]) # BELT
AiL = np.array([[0.,-1.,0.,0.],[-1.,0.,0.,0.],[0.,0.,-1.,0.],[0.,0.,0.,1.]]) # WHEEL

lead_dist = 60
blade = add_lead_points(blade, lead_dist)
cx_frenes = get_cx_frenes(blade, 2, 3)

cx_path = get_cx_path(cx_frenes, AiL)
rounded_cx_path = np.round(np.array(cx_path), 3)

generate_cx_dat("convex_wheel", "A", rounded_cx_path)
generate_cx_src("convex_wheel", "A", rounded_cx_path)


'''
### BACKUP ###

# T -> 0
oT = np.array([1009.15, -16.49, 623.81])
xT = np.array([996.14, 1010.89, 1010.89, 1023.99, 1014.15, 1014.15, 1004.89, 1004.89, 1009.15])
yT = np.array([-16.14, -29.24, 0.92, -16.14, -10.54, -22.95, -22.21, -10.51, -16.49])
zT = np.array([625.57, 623.52, 623.48, 622.35, 623.61, 622.86, 624.73, 624.40, 623.81])
[BELT_FRAME, AT0] = get_belt_frame(oT, xT, yT, zT)
#print("BELT_FRAME:\n", "X: ", BELT_FRAME[0], "| Y: ", BELT_FRAME[1], "| Z: ", BELT_FRAME[2], "| A: ", BELT_FRAME[3], "| B: ", BELT_FRAME[4], "| C: ", BELT_FRAME[5])

class Blade:
    def __init__(self, data: dict):
        self._convex = []
        self._npc = [] # number of points on the convex/concave
        for key, value in data.items():
            x_cx = value["x_cx"]; y_cx = value["y_cx"]; z = value["z"]
            self._npc.append(len(x_cx))
            self._convex.append([])
            for x,y in zip(x_cx, y_cx):
                self._convex[-1].append(np.array([x,y,z]))
    @property
    def convex(self):
        return self._convex
    @property
    def npc(self):
        return self._npc

def get_lead_point(frene1: Frene, frene2: Frene, dist: float) -> Frene:
    p1 = frene1.transf[0:3,3].transpose()
    p2 = frene2.transf[0:3,3].transpose()
    v21 = p1[0:3] - p2[0:3]
    v13 = dist/np.linalg.norm(v21)*v21
    p3 = p1[0:3] + v13
    return Frene(frene1.t, frene1.b, frene1.n, p3)

def get_airstripe_frenes(i: int, blade: Blade) -> list:
    frenes = []
    for ptnum in range(0,blade.npc[i]): #0, blade.npc[i]
        if(ptnum == 0):
            p0 = blade.convex[i][0]
            u1 = blade.convex[i][1]
            u2 = blade.convex[i][2]
        elif(ptnum == blade.npc[i]-1):
            p0 = blade.convex[i][blade.npc[i]-1]
            u1 = blade.convex[i][blade.npc[i]-2]
            u2 = blade.convex[i][blade.npc[i]-3]
        else:
            p0 = blade.convex[i][ptnum]
            u1 = blade.convex[i][ptnum-1]
            u2 = blade.convex[i][ptnum+1]   
        v1 = blade.convex[i-1][ptnum]
        frenes.append(get_frene(p0, u1, u2, v1))        
    return frenes

def get_blade_frames(path):
    with open(path, 'r') as file:
        blade = Blade(json.load(file))
    return get_airstripe_frenes(1, blade)
'''