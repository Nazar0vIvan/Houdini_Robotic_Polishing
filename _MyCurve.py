from hou.session import *
import hou

assetNode = kwargs['type'].instances()[0] # asset node object
drawCurveNode = assetNode.children()[0]

parmNodePointsInput = assetNode.parms()[0]
parmNodePoints = assetNode.parms()[1]
if(parmNodePointsInput.eval() == 0):
    parmNodePoints.disable(True)
# ----------------------------------------------------------------------------------------
def houPointsToArray(houPoints):
    nodePoints = []
    for houPoint in houPoints:
        nodePoints.append([houPoint.position().x(),
                           houPoint.position().y(),
                           houPoint.position().z()])
    return np.array(nodePoints)
# ----------------------------------------------------------------------------------------        
def parsePointsCoordinates():
    nodePoints = []
    pointsStr = parmNodePoints.eval().split(' ')
    for pointStr in pointsStr:
        nodePoints.append(map(float, pointStr.split(',')))
    return nodePoints
# ----------------------------------------------------------------------------------------
def changeNodePointsInputAndRedrawCurve():
    input_geo = assetNode.inputs()[0].geometry() # input node geometry
    parmNodePointsInput = assetNode.parms()[0]
    parmNodePoints = assetNode.parms()[1]
    if (parmNodePointsInput.eval() == 0): # By Input
        parmNodePoints.disable(True)
        if (len(input_geo.points()) != 0):
            curve.setNodePoints(houPointsToArray(input_geo.points()))
            curve.interpolate()
        else:
            hou.ui.displayMessage("There are no node points in input geometry")
    elif(parmNodePointsInput.eval() == 1): # By Hand
        parmNodePoints.disable(False)
        if (parmNodePoints.eval() != ""):
            curve.setNodePoints(parsePointsCoordinates())
            curve.interpolate()
        else:
            hou.ui.displayMessage("Please, enter node points")
    drawCurveNode.cook(force=True)
# ----------------------------------------------------------------------------------------    
def changeInterpolationStrategyAndRedrawCurve(): 
    input_geo = assetNode.inputs()[0].geometry() # input node geometry
    parmInterpolationStrategy = assetNode.parms()[2]  
    if (parmInterpolationStrategy.eval() == 0):  # LIN
        curve.setInterpolationStrategy(LinearInterpolation())
    elif(parmInterpolationStrategy.eval() == 1): # CIRC
        curve.setInterpolationStrategy(CircInterpolation())
    elif(parmInterpolationStrategy.eval() == 2): # CUBIC
        curve.setInterpolationStrategy(CubicInterpolation())
    elif(parmInterpolationStrategy.eval() == 3): # QUADRATIC
        curve.setInterpolationStrategy(QuadraticInterpolation())
    elif(parmInterpolationStrategy.eval() == 4): # POLYNOMIAL
        curve.setInterpolationStrategy(PolyInterpolation())
    curve.interpolate()
    drawCurveNode.cook(force=True)
# ----------------------------------------------------------------------------------------    
def changeStep():
    parmStep = assetNode.parms()[3]
    curve.setStep(parmStep.evalAsFloat())
    curve.interpolate()
    drawCurveNode.cook(force=True)
# ----------------------------------------------------------------------------------------
def getPoints():
    return curve.calculateAndGetPoints()
# ----------------------------------------------------------------------------------------
    print("input changed")#On Updated works when you jump to this node from other node

#assetNode = kwargs['node']
#asset = assetNode.hdaModule()

#asset.changeNodePointsInputAndRedrawCurve()




from hou.session import *

node = hou.pwd()
geo = node.geometry()

geo.deletePoints(geo.points())

curveAsset = hou.node("../../curve/").hdaModule()
points = curveAsset.getPoints()

if(len(points) !=0):
    for point in points:
        p = geo.createPoint()
        p.setPosition((point[0],point[1],point[2]))
    for i in range(len(geo.points())-1):
        poly = geo.createPolygon()
        poly.addVertex(geo.points()[i])
        poly.addVertex(geo.points()[i+1])        
        
        
