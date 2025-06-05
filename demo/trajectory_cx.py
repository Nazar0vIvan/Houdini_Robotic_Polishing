import hou
import numpy as np
from math import *
import sympy as sym
from linalg import *
from pathplanner import *
from tools import *
from diffgeom import *
from blade import * 

# -----
def draw_frene(geo, p, frene, scale=10):
    pt = geo.createPoint()
    pt.setPosition((p[0], p[1], p[2]))
    
    tpt = geo.createPoint()
    tpt.setPosition((float(p[0]+scale*frene[0][0]), float(p[1]+scale*frene[1][0]), float(p[2]+scale*frene[2][0])))
    tpl = geo.createPolygon()
    tpl.addVertex(pt)
    tpl.addVertex(tpt)
    
    bpt = geo.createPoint()
    bpt.setPosition((float(p[0]+scale*frene[0][1]), float(p[1]+scale*frene[1][1]), float(p[2]+scale*frene[2][1])))
    bpl = geo.createPolygon()
    bpl.addVertex(pt)
    bpl.addVertex(bpt)
    
    npt = geo.createPoint()
    npt.setPosition((float(p[0]+scale*frene[0][2]), float(p[1]+scale*frene[1][2]), float(p[2]+scale*frene[2][2])))
    npl = geo.createPolygon()
    npl.addVertex(pt)
    npl.addVertex(npt)
# -----

# i - frame, attached to the point on the airfoil surface
# S - frame, attached to the point on the roll surface
# B - blade frame
# C - roll frame (center of the roll)

# Tr_i_j - from i to j
# MAIN FOLMULA: Tr_C_0 = Tr_i_0 * Tr_S_i * Tr_C_S

np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)

def create_animation():
    hou.playbar.addEventCallback(robotMoveEvent)

def robotMoveEvent(event_type, frame):
    tool_obj = hou.node("obj/1FF1_150x16x32")
    setToolPosition(tool_obj, path_planner.path[frame])

def setToolPosition(node, position):
    node.parm("tx").set(position[0])
    node.parm("ty").set(position[1])
    node.parm("tz").set(position[2])
    node.parm("rx").set(position[3])
    node.parm("ry").set(position[4])
    node.parm("rz").set(position[5])

def get_frene(p0, u1, u2, v1):
    a0,a1,a2,t = sym.symbols('a0 a1 a2 t')
    xu = t; yu = a0*t**2+a1*t+a2; zu = u1[2];             
    vecfun_u = VectorFunction(xu, yu, zu)
    coeffs_u = poly(u1[0], p0[0], u2[0], u1[1], p0[1], u2[1]) 
    tanu = vecfun_u.tangent_val(p0[0], coeffs_u)
    if(tanu[0] > 0): tanu = -1*tanu 
    
    tanv = v1 - p0
    normalized_tanv = tanv/np.linalg.norm(tanv)

    norm = np.cross(tanu, normalized_tanv)
    normalized_norm = norm/np.linalg.norm(norm)
    binorm = np.cross(normalized_norm, tanu)
    
    return Frene(tanu, binorm, norm, p0)

# ------------------- path functions ------------------- # 

def get_lead_frene(frene1: Frene, frene2: Frene, dist: float) -> Frene:
    p1 = frene1.transf[0:3,3].transpose()
    p2 = frene2.transf[0:3,3].transpose()
    v21 = p1[0:3] - p2[0:3]
    v13 = dist/np.linalg.norm(v21)*v21
    p3 = p1[0:3] + v13
    return Frene(frene1.t, frene1.b, frene1.n, p3)

def get_airstripe_path(i: int, blade: Blade, tool: Tool, lead_in_dist: float = 0.0, lead_out_dist: float = 0.0) -> list:
    frenes = []
    for ptnum in range(0, blade.npc): #0,blade.npc
        if(ptnum == 0):
            p0 = blade.convex[i][0]
            u1 = blade.convex[i][1]
            u2 = blade.convex[i][2]
        elif(ptnum == blade.npc-1):
            p0 = blade.convex[i][blade.npc-1]
            u1 = blade.convex[i][blade.npc-2]
            u2 = blade.convex[i][blade.npc-3]
        else:
            p0 = blade.convex[i][ptnum]
            u1 = blade.convex[i][ptnum-1]
            u2 = blade.convex[i][ptnum+1]   
        v1 = blade.convex[i-1][ptnum]
        frenes.append(get_frene(p0,u1,u2,v1))        
    lead_in_frene = get_lead_frene(frenes[0], frenes[1], lead_in_dist)
    lead_out_frene = get_lead_frene(frenes[-1], frenes[-2], lead_out_dist)
    frenes.insert(0, lead_in_frene)
    frenes.append(lead_out_frene)
    return frenes2path(frenes, tool)

def frenes2path(frenes: list, tool: Tool) -> dict:
    tool_path = []; transforms = []
    for frene in frenes:
        tool.transf_orig2base = np.dot(frene.transf, tool.transf_orig2blank)
        #transforms.append(tool_transf) # if transformation matrices are needed
        euler = rot2euler(tool.transf_orig2base[0:3,0:3], True)
        tool_path.append([*tool.transf_orig2base[0:3,3], euler["C1"], euler["B1"], euler["A1"]])
    return tool_path

def attrpoints2dict(x: list, y: list, z: list, npf: int, npt: int):
    result = {}
    for i in range(npf):
        result[i] = []
        for j in range(npt):
            result[i].append(np.array([x[i*npt+j], y[i*npt+j], z[i*npt+j]]))
    return result 
    
def get_airstripes_count(cut_depth: float, tool: Tool, blade: Blade):
    cut_width = 2*sqrt(cut_depth*(2*tool.B/2 - cut_depth))
    return 2*ceil(blade.airheight/cut_width)

def save_path(path: list):
    with open('path.txt', 'w') as file:
        file.writelines("[" + ", ".join("{:6.3f}".format(el) for el in pos) + "]\n" for pos in path)
    file.close()

def path2cls(feed, speed, r, path = []):
    # feed - m/min
    # speed - m/s
    # r - mm
    with open("tool_trajectory.cls", 'w') as file:
        file.write(f"FEDRAT/MMPM,{feed*1000:.0f}\n")
        file.write(f"SPINDL/RPM,{speed*1000.*60./(2.*np.pi*r):.0f},CLW\n")
        for p in path:
            str_p = ','.join(f'{el:.4f}' for el in p)
            file.write(f"GOTO/{str_p}\n")

# ------------------- START ------------------- #    

# --------------- initial data ---------------- #
feed = 1.0 # m/min
speed = 20 # m/s [n = v/(2*pi*R)]
cut_depth = 0.03 # mm
lead_in_dist = 20 # mm
lead_out_dist = 20 # mm
stripes_to_animate = 2
# --------------------------------------------- #

hou.setFps(24)

wheel_air = W1FF1(150.0, 16.0)
wheel_air.transf_surf2blank = np.array([[0.,0.,-1.,0.], [0.,-1.,0.,0.], [-1.,0.,0.,0.], [0.,0.,0.,1.]])

blade142 = hou.node("/obj/blade142/blade142")

npf = blade142.geometry().intAttribValue("pf_count")
npc = blade142.parm("pt_cxcv_count").eval()
npe = blade142.parm("pt_edges_count").eval()
npr = blade142.parm("pt_radius_count").eval()

x_cx = list(blade142.geometry().floatListAttribValue("x_cx"))
y_cx = list(blade142.geometry().floatListAttribValue("y_cx"))
z_cx = list(blade142.geometry().floatListAttribValue("z_cx"))

convex = attrpoints2dict(x_cx, y_cx, z_cx, npf, npc)

blade = Blade(convex, {}, npf, npc, npe, npr)

stipes_count = get_airstripes_count(cut_depth, wheel_air, blade)
blade142.parm("pf_count").set(stipes_count)

path_planner = PathPlanner(24)
cut_path = []
for pfnum in range(1, stripes_to_animate): #(1, blade.npf)
    airstripe_basepath = get_airstripe_path(pfnum, blade, wheel_air, lead_in_dist, lead_out_dist)
    stripe_time = blade.pf_cx_length(pfnum)/feed
    stripe_extra_points_count = floor(stripe_time*hou.fps())
    extra_points_count = floor((stripe_extra_points_count - blade.npc)/(blade.npc-1))
    for i in range(1,len(airstripe_basepath)-2):
        cur_pos = airstripe_basepath[i]
        next_pos = airstripe_basepath[i+1]
        cut_path_segment = np.linspace(cur_pos, next_pos, extra_points_count).tolist()
        cut_path_segment.pop()
        cut_path.extend(cut_path_segment) 
    path_planner.add_lin(airstripe_basepath[0], airstripe_basepath[1], feed)
    path_planner.add_path(cut_path)
    path_planner.add_lin(airstripe_basepath[-2], airstripe_basepath[-1], feed)
    cut_path = []

hou.playbar.setFrameRange(0, len(path_planner.path)-1)

create_animation()
#hou.playbar.clearEventCallbacks()

path2cls(feed, speed, wheel_air.D/2, path_planner.path)
