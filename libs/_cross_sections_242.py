import blade
import ruledsurf

# -------------------------- TRAJECTORY -------------------------- #

# i - frame, attached to the point on the airfoil surface
# S - frame, attached to the point on the roll surface
# B - blade frame
# C - roll frame (center of the roll)

np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)

def create_animation():
    hou.playbar.addEventCallback(robotMoveEvent)

def robotMoveEvent(event_type, frame):
    roll_obj = hou.node("obj/roll")
    setToolPosition(roll_obj, trajectory[0][frame])

def setToolPosition(node, position):
    node.parm("tx").set(position[0])
    node.parm("ty").set(position[1])
    node.parm("tz").set(position[2])
    node.parm("rx").set(position[3])
    node.parm("ry").set(position[4])
    node.parm("rz").set(position[5])

def get_ruled_surf(blade_data: dict, i: int, j: int, nu: int, nv: int) -> RuledSurf:
    a1 = blade.geo[i][j];   a2 = blade.geo[i][j+1]
    b1 = blade.geo[i+1][j]; b2 = blade.geo[i+1][j+1]
    return RuledSurf(a1,a2,b1,b2,nu,nv)
    
def get_ruled_surface_frenes(ruled_surf: RuledSurf) -> list:   
    step_v = 1./ruled_surf.nv
    step_u = 1./ruled_surf.nu
    v = 0
    frenes = dict.fromkeys(np.arange(ruled_surf.nv), [None]*ruled_surf.nu)
    for iv in range(ruled_surf.nv):
        u = 0.
        for iu in range(ruled_surf.nu):
            pt = ruled_surf.val(u,v)
            du = ruled_surf.du(v)
            dv = -ruled_surf.dv(u)
            binorm = dv
            normal = np.cross(du,dv)
            tangent = np.cross(binorm, normal)
            frenes[iv][iu] = np.column_stack((
              np.append(tangent,0.),
              np.append(binorm,0.),
              np.append(normal,0.),
              np.append(pt,1.))
            )
            u = u + step_u
        v = v + step_v
    return frenes
    
def update_trajectory(trajectory: dict, frenes: dict, Tr: np.array, i: int) -> dict:
    nu = len(frenes[0])
    for iv in frenes.keys():
        for iu in range(nu):
            Tr_C_0 = np.dot(frenes[iv][iu], Tr)
            position = Tr_C_0[0:3,3]
            euler = rot2euler(Tr_C_0[0:3,0:3], True)
            trajectory[iv][i*nu + iu] = np.array([
                Tr_C_0[0,3],
                Tr_C_0[1,3],
                Tr_C_0[2,3],
                euler["C1"],
                euler["B1"],
                euler["A1"]]
            )
    return trajectory
     
def get_trajectory(blade: dict, velocity: float, r: float, h: float) -> list:  
  
    nv = floor(blade.airfoil_height/h)+1
    
    time_cx = blade.cx_length/velocity
    frames_count_cx = floor(time_cx*hou.fps())
    nu = floor(frames_count_cx/(blade.npc-1))
    
    Tr_S_i = np.array([[0.,0.,1.,0.], [0.,1.,0.,0.], [-1.,0.,0.,0.], [0.,0.,0.,1.]])
    Tr_C_S = np.array([[1.,0.,0.,-r], [0.,1.,0.,h/2.], [ 0.,0.,1.,0.], [0.,0.,0.,1.]])
    Tr_C_i = np.dot(Tr_S_i, Tr_C_S)
    
    trajectory = dict.fromkeys(np.arange(nv), [None]*nu*nv*(blade.npc-1))
    for i in range(1): #blade.npc-1
        ruled_surf = get_ruled_surf(blade, 0, i, nu, nv)
        frenes = get_ruled_surface_frenes(ruled_surf)
        trajectory = update_trajectory(trajectory, frenes, Tr_C_i, i)
    return trajectory, nv*frames_count_cx

blade142 = hou.node("/obj/blade142/blade142")
npf = blade142.geometry().intAttribValue("pf_count")
npc = blade142.parm("pt_cxcv_count").eval()
npe = blade142.parm("pt_edges_count").eval()
npt = 2*npc+2*npe

x = list(blade142.geometry().floatListAttribValue("x"))
y = list(blade142.geometry().floatListAttribValue("y"))
z = list(blade142.geometry().floatListAttribValue("z"))

blade = Blade(x, y, z, npf, npc, npe)

roll = hou.node("/obj/roll/roll")
r = roll.parmTuple("rad").eval()[0]
h = roll.parm("height").eval()
        
hou.setFps(24)

trajectory, frame_count = get_trajectory(blade, 1, r, h)

hou.playbar.setFrameRange(0, frame_count)
create_animation()
