import numpy as np
import json
from linalg import *
from diffgeom import *
from math import radians

np.set_printoptions(suppress=True)
np.set_printoptions(precision=3)

def print_float_lists(list_of_lists):
    print("[")
    for inner_list in list_of_lists:
        formatted = [float(f"{num:.3f}") for num in inner_list]
        print(f"    {formatted},")
    print("]")

### BELT ###

def get_belt_frame(orig: np.array, x: np.array, y: np.array, z: np.array):
    # T - TOOL: belt - tau, beta, nu

    [A, B, C, D, AA, BB, DD] = points2plane(x, y, z)

    # belt normal - z
    nu = np.array([A, B, C])

    # belt tau - x
    tau_xend = orig[0]-40
    tau_yend = orig[1]
    tau_zend = AA*(tau_xend) + BB*(tau_yend) + DD
    tau_end = np.array([tau_xend, tau_yend, tau_zend])
    tau = tau_end - orig
    mtau = np.linalg.norm(tau)
    tau=tau/mtau

    # belt beta - y
    beta = np.cross(nu,tau)

    # euler
    tau = tau.reshape((-1,1))
    beta = beta.reshape((-1,1))
    nu = nu.reshape((-1,1))
    RT0 = np.hstack((tau, np.hstack((beta,nu)))) # rotation matrix from belt frame to robot base frame

    A = np.degrees(np.arctan2(RT0[1,0],RT0[0,0]))
    B = np.degrees(np.arcsin(-RT0[2,0]))
    C = np.degrees(np.arctan2(RT0[2,1],RT0[2,2]))
    frame = np.array([orig[0], orig[1], orig[2], A, B, C])

    # matrix AT0
    AT0 = RT0
    AT0 = np.append(RT0, orig.reshape((-1,1)), axis=1)
    AT0 = np.vstack((AT0, np.array([0,0,0,1])))

    return [frame, AT0]

### BLADE ### 

def get_lead_point(pt1: list, pt2: list, dist: float) -> np.ndarray:
    #pt1 - the last point of cx/cv
    #pt2 - the point before the last of cx/cv
    pt1 = np.array(pt1)
    pt2 = np.array(pt2)
    v21 = pt1 - pt2
    v13 = dist/np.linalg.norm(v21)*v21
    return list(np.round(pt1 + v13, 3))

def get_frene(p0: np.array, u1: np.array, u2: np.array, v1: np.array) -> Frene:
    a0,a1,a2,t = sym.symbols('a0 a1 a2 t')
    xu = t; yu = a0*t**2+a1*t+a2; zu = u1[2];             
    vecfun_u = VectorFunction(xu, yu, zu)
    coeffs_u = poly(u1[0], p0[0], u2[0], u1[1], p0[1], u2[1]) 
    tanu = vecfun_u.tangent_val(p0[0], coeffs_u)
    if(tanu[0] < 0): tanu = -1*tanu
    tanv = v1 - p0
    normalized_tanv = tanv/np.linalg.norm(tanv)
    norm = np.cross(tanu, normalized_tanv)
    normalized_norm = norm/np.linalg.norm(norm)
    binorm = np.cross(normalized_norm, tanu)
    return Frene(tanu, binorm, norm, p0)

# CX - BEGIN
def get_cx_frenes(blade: dict, airstripes_count: int) -> list:
    frenes = []
    for i in range(airstripes_count):
        frenes.append(get_cx_airstripe_frenes(blade, i))
    return frenes

def get_cx_airstripe_frenes(blade: dict, pfnum: int) -> list:
    pf_frenes = [Frene()]
    npt = len(blade[pfnum]["cx"])
    for ptnum in range(1, npt-1):
        p0 = np.array(blade[pfnum]["cx"][ptnum])
        u1 = np.array(blade[pfnum]["cx"][ptnum-1])
        u2 = np.array(blade[pfnum]["cx"][ptnum+1])   
        v1 = np.array(blade[pfnum+1]["cx"][ptnum])
        pf_frenes.append(get_frene(p0, u1, u2, v1))
    pf_frenes[0] = Frene(pf_frenes[1].t, pf_frenes[1].b, pf_frenes[1].n, np.array(blade[pfnum]["cx"][0]))
    pf_frenes.append(Frene(pf_frenes[-1].t, pf_frenes[-1].b, pf_frenes[-1].n, np.array(blade[pfnum]["cx"][-1])))
    return pf_frenes if pfnum % 2 == 0 else pf_frenes[::-1]

def generate_cx_dat(program_name: str, pf_prefix: str, path: list) -> None:
    n = len(path)
    m = len(path[0])
    with open(f"{program_name}.dat", 'w') as file:
        file.write(f"DEFDAT {program_name}\n")
        file.write("EXTERNAL DECLARATION")
        for i in range(n):
            file.write(f"\n\n; ### {pf_prefix}{i+1} ###\n\n")
            file.write(f"DECL POS {pf_prefix}{i+1}[{m}]\n\n")
            for j, pos in enumerate(path[i]):
                file.write(f"{pf_prefix}{i+1}[{j+1}] = {{X {pos[0]}, Y {pos[1]}, Z {pos[2]}, A {pos[3]}, B {pos[4]}, C {pos[5]}}}\n")
            
def generate_cx_src(program_name: str, pf_prefix: str, path: list) -> None:
    n = len(path)
    m = len(path[0])
    with open(f"{program_name}.src", 'w') as file:
        file.write(f"DEF {program_name}()\n")
        file.write("INI\n\n")
        file.write("PTP HOME Vel=100% DEFAULT\n\n")
        file.write("$BASE = BASE_DATA[6]\n")
        file.write("$TOOL = TOOL_DATA[2]\n\n")
        file.write("PTP {A1 0, A2 -90, A3 90, A4 0, A5 0, A6 0}")
        for i in range(n):
            file.write(f"\n\n; ### {pf_prefix}{i+1} ###\n")
            file.write("SPLINE WITH $VEL= {CP 0.001, ORI1 0.001, ORI2 0.001}\n")
            for j, pos in enumerate(path[i]):
                file.write(f"  SPL {pf_prefix}{i+1}[{j+1}]\n")
            file.write("ENDSPLINE")

        
# CX - END

### BEGIN ###

with open("99.01.25.242.json", 'r') as file:
    blade = json.load(file)

# i -> B
cx_frenes = []
'''
[
    [Frene_11, Frene_12 ...]
    [Frene_21, Frene_22 ...]
    ...
    [Frene_n1, Frene_n2 ...]
]
'''


for i in range(2):
    blade[i]["cx"][0] = get_lead_point(blade[i]["cx"][1], blade[i]["cx"][2], 5)
    blade[i]["cx"][-1] = get_lead_point(blade[i]["cx"][-2], blade[i]["cx"][-3], 5)

cx_frenes = get_cx_frenes(blade, 2)

# B -> F
ABF = np.dot(translationMatrix([0.011, 0.047, 153.319]), rotationMatrix4x4(radians(-49), "z"))
# i -> T
AiT = np.array([[0,-1,0,0],[-1,0,0,0],[0,0,-1,0],[0,0,0,1]])

ABTs = []
path = []
#i = 0; j = 0
# ABF_inv = np.linalg.inv(ABF)
for pf_frene in cx_frenes:
    pf_path = []
    #j = 0
    for frene in pf_frene:
        ABT = np.dot(AiT, np.linalg.inv(frene.transf))
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
        
rounded_path = np.round(np.array(path),3)

generate_cx_dat("convex_wheel", "A", rounded_path)
generate_cx_src("convex_wheel", "A", rounded_path)

'''
i = 0
for pf in rounded_path:
    print(f"profile #{i}")
    for p in pf:
        print(    p)
    i = i + 1
'''


'''
def printFloatList(data, precision):
    for i in range(len(data)):
        print [f"%0.{precision}f" % i for i in a]
'''

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