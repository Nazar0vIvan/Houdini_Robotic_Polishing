import socket
import json
#from kinematics import *
from hou.session import *

trajectory = []

def slotTest(kwargs,assetNode):
    print(assetNode.parm("shoulder").eval())
# ---------------------------------------------------------------------------------------
def updateAF0s(kwargs,assetNode):

    ABF_0 = assetNode.parmTuple("ABF_0").eval()
    ABF_1 = assetNode.parmTuple("ABF_1").eval()
    ABF_2 = assetNode.parmTuple("ABF_2").eval()
    ABF_3 = assetNode.parmTuple("ABF_3").eval()
    
    AT0_0 = assetNode.parmTuple("AT0_0").eval()
    AT0_1 = assetNode.parmTuple("AT0_1").eval()
    AT0_2 = assetNode.parmTuple("AT0_2").eval()
    AT0_3 = assetNode.parmTuple("AT0_3").eval()
    
    AiT_0 = assetNode.parmTuple("AiT_0").eval()
    AiT_1 = assetNode.parmTuple("AiT_1").eval()
    AiT_2 = assetNode.parmTuple("AiT_2").eval()
    AiT_3 = assetNode.parmTuple("AiT_3").eval()
    
    ABF = np.array([ABF_0,ABF_1,ABF_2,ABF_3])
    AT0 = np.array([AT0_0,AT0_1,AT0_2,AT0_3])
    AiT = np.array([AiT_0,AiT_1,AiT_2,AiT_3])
    
    blade = hou.node(assetNode.parm("trajectorySOP").eval()) # get blade node
    AiBs_c = blade.parm("AiBs_c").eval()
    AiBs_c = json.loads(AiBs_c)

    AF0s = {"values": []}
    convex = AiBs_c.get("convex")
    for AiBs_c in convex.values():
        for AiB_c in AiBs_c:
            AiB_c = np.array(AiB_c)
            AF0 = np.dot(AT0,np.dot(AiT,np.dot(np.linalg.inv(AiB_c),np.linalg.inv(ABF))))
            AF0s.get("values").append(AF0.tolist())  
    assetNode.parm("AF0s").set(json.dumps(AF0s))

### ----------------------------------- PARAMETERS ---------------------------------- ###
def updateAT0(kwargs,assetNode):
    toolNodePath = assetNode.parm("toolNodePath").eval()
    if(not toolNodePath): return
    grinder = hou.node(toolNodePath)
    assetNode.parmTuple("AT0_0").set(grinder.parmTuple("ATS0_0").eval())
    assetNode.parmTuple("AT0_1").set(grinder.parmTuple("ATS0_1").eval())
    assetNode.parmTuple("AT0_2").set(grinder.parmTuple("ATS0_2").eval())
    updateAF0s(kwargs,assetNode)
# ---------------------------------------------------------------------------------------

### ------------------------------------ CONTROLS ----------------------------------- ###
def slotTrajectoryJSONChanged(kwargs,assetNode):
    try:
        f = open(kwargs["parm"].eval())
        data = json.load(f) # trajectory json
    except json.decoder.JSONDecodeError:
        hou.ui.displayMessage("Decoding JSON has failed")
# ---------------------------------------------------------------------------------------
def slotTrajectorySOPChanged(kwargs,assetNode):
    if(not kwargs["parm"].eval()): return
    updateAF0s(kwargs,assetNode)
# ---------------------------------------------------------------------------------------
def slotStartSim(kwargs,assetNode):
    AF0s = json.loads(assetNode.parm("AF0s").eval())
    if (not AF0s): hou.ui.displayMessage("There are no points in the trajectory")
    
    global trajectory 
    trajectory = []
    AF0s = AF0s.get('values') 
    for AF0 in AF0s: 
        trajectory.append(np.array(AF0))
        #x = AF0[0][3]
        #y = AF0[1][3]
        #z = AF0[2][3]
        #euler = rotationMatrixToEulerAngles(AF0)
        #a = euler.get('A1')
        #b = euler.get('B1')
        #c = euler.get('C1')
        # trajectory.append(AF0)
        # trajectory.append([x,y,z,a,b,c]) # !global!
    
    npoints = len(trajectory)
    T = (npoints-1)*4.0/1000.0 # total time to complete the path
    hou.setFps(ceil(npoints/T))
    hou.playbar.setFrameRange(0,npoints-1)
    hou.playbar.addEventCallback(robotMoveEvent)
    # hou.playbar.play()
        
'''
    trajectoryPointsFromGeo = hou.node(parmTrajectoryPath.eval()).geometry().points()
        
    if (len(trajectoryPointsFromGeo) == 0):
        hou.ui.displayMessage("There are no points in the trajectory")
    else:            
        trajectoryObjNode = hou.node(parmTrajectoryPath.eval()).parent() 
        ATW_hou = trajectoryObjNode.worldTransform().transposed().asTupleOfTuples()
        ATW = np.asarray(ATW_hou)
        ATW[np.absolute(ATW)<=0.0001]=0
        
        A0W_hou = nodeBase.worldTransform().transposed().asTupleOfTuples()
        A0W = np.asarray(A0W_hou)
        A0W[np.absolute(A0W)<=0.0001]=0
        
        A0i = np.dot(np.linalg.inv(A0W),ATW)
        global trajectoryPoints
        trajectoryPoints = []
        for point in trajectoryPointsFromGeo:
            p = np.dot(A0i,np.array([[point.position().x()],[point.position().y()],[point.position().z()],[1]]))
            trajectoryPoints.append(p)      
        npoints = len(trajectoryPointsFromGeo)
        T = (npoints-1)*4.0/1000.0 # total time to compelete the path
        hou.setFps(ceil(npoints/T))
        hou.playbar.setFrameRange(0,npoints-1)
        hou.playbar.addEventCallback(robotMoveEvent)
#            hou.playbar.play()
        print("you can start simulation")
'''
# ----------------------------------------------------------------------------------------
def slotStopSim(kwargs,assetNode):
    hou.playbar.stop()
    hou.playbar.clearEventCallbacks()
    print("simulation is stopped")
# ----------------------------------------------------------------------------------------
def ik(AF0, type = "manual"):    
    assetNode = kwargs['type'].instances()[0]

    if(type == "manual"):
        x = assetNode.parm("X").eval(); y = assetNode.parm("Y").eval(); z = assetNode.parm("Z").eval()
        a = assetNode.parm("A").eval(); b = assetNode.parm("B").eval(); c = assetNode.parm("C").eval()

    shoulder = assetNode.parm("shoulder").eval()
    elbow = assetNode.parm("elbow").eval()

    solution = solveIK(AF0, shoulder, elbow)
    
    qt = solution.get("qt")
    qr = solution.get("qr")
    shoulder = solution.get("shoulder")
    elbow = solution.get("elbow")
    
    #print(qt)
    #print([qr[0][0],qr[1][0],qr[2][0]])
    #print('---')
    
    assetNode.parm("shoulder").set(shoulder)
    assetNode.parm("elbow").set(elbow)
    assetNode.parm("q1").set(qt[0])
    assetNode.parm("q2").set(qt[1])
    assetNode.parm("q3").set(qt[2])
    assetNode.parm("q4").set(qr[0][0])
    assetNode.parm("q5").set(qr[1][0])
    assetNode.parm("q6").set(qr[2][0])
    #assetNode.parmTuple("H1").set(AF0[0,0:4])
    #assetNode.parmTuple("H2").set(AF0[1,0:4])
    #assetNode.parmTuple("H3").set(AF0[2,0:4])
    #assetNode.parmTuple("H4").set([0,0,0,1])
# ---------------------------------------------------------------------------------------
def dk():
    # dk solution is unusual, dk is needed when we change q vector and 
    # while it happens we obtain flange -> base transformation matrix A6B
    # from flange node and then derived X,Y,Z,A,B,C from it  
    assetNode = kwargs['type'].instances()[0]
    
    BASE = assetNode.node("0_BASE")
    FLANGE = assetNode.node("FLANGE")
    q2 = assetNode.parm("q2").eval()
    q3 = assetNode.parm("q3").eval()
    
    solution = solveDK(BASE,FLANGE,q2,q3) # return dict "TCP", "A6B", "shoulder", "elbow"
    tcp = solution.get("TCP")
    AF0 = solution.get("A6B")
    
    assetNode.parm("X").set(tcp.X)
    assetNode.parm("Y").set(tcp.Y)
    assetNode.parm("Z").set(tcp.Z)
    assetNode.parm("A").set(tcp.A)
    assetNode.parm("B").set(tcp.B)
    assetNode.parm("C").set(tcp.C)
    assetNode.parm("shoulder").set(solution.get("shoulder"))
    assetNode.parm("elbow").set(solution.get("elbow"))
    assetNode.parmTuple("H1").set(AF0[0,0:4])
    assetNode.parmTuple("H2").set(AF0[1,0:4])
    assetNode.parmTuple("H3").set(AF0[2,0:4])
    assetNode.parmTuple("H4").set([0,0,0,1])
# ---------------------------------------------------------------------------------------
def robotMoveEvent(event_type,frame):
    # for each frame ik will be solved
    #x = trajectory[frame][0]
    #y = trajectory[frame][1]
    #z = trajectory[frame][2]
    #a = trajectory[frame][3]
    #b = trajectory[frame][4]
    #c = trajectory[frame][5]
    ik(trajectory[frame],"not manual")
# ----------------------------------------------------------------------------------------
def slotOpenPort(kwargs,assetNode):
    print("open port")
# ----------------------------------------------------------------------------------------
def slotClosePort(kwargs,assetNode):
    print("close port") 
# ----------------------------------------------------------------------------------------
def home():
    assetNode = kwargs['type'].instances()[0]
    assetNode.parm("q1").set(0.0)
    assetNode.parm("q2").set(-90.0)
    assetNode.parm("q3").set(90.0)
    assetNode.parm("q4").set(0.0)
    assetNode.parm("q5").set(0.0)
    assetNode.parm("q6").set(0.0)  
    dk()