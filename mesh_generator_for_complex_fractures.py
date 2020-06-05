import math
import random
import datetime
import os

###################################################
###             USER INPUTS                     ###
###################################################
isMesh3D = False
dynamicCrack = False
staticCrack = not dynamicCrack
glueTol = 0.05
tolResBB = 0.1
useQuadElementsRes = 0
useQuadElementsFrac = 0
maxGridSize = 5
resOuterGridSize = maxGridSize
fracGridMinSize = 1.0
fracGridMaxSize = fracGridMinSize
resGridMinSize = fracGridMinSize
minGridSize = fracGridMinSize
SetFineness3D = 0     # 0 for very coarse; 1 for very coarse; 2 for moderate; 3 for fine; 4 for very fine; 5 for custom;
SetGrowthRate3D = 0.2   # used only when SetFineness3D = 5
SetNbSegPerEdge3D = 2   # used only when SetFineness3D = 5
SetNbSegPerRadius3D = 2 # used only when SetFineness3D = 5

# for 2D mesh
SetFineness2D = 0     # 0 for very coarse; 1 for very coarse; 2 for moderate; 3 for fine; 4 for very fine; 5 for custom;
SetGrowthRate2D = 0.2     # used only when SetFineness2D = 5
SetNbSegPerEdge2D = 2     # used only when SetFineness2D = 5
SetNbSegPerRadius2D = 2   # used only when SetFineness2D = 5  
zGridSize = 5
timesExtrusion = 1

# reservoir dimensions
fracLengthScaling = 0.5      # max frac length relative to reservoir dimensions
fracCenterExtent = 0.4       # location of fracture centre w.r.t. reservoir dimension
xNegExtra = 0
xPosExtra = xNegExtra
yNegExtra = 0
yPosExtra = yNegExtra
zNegExtra = 0
zPosExtra = zNegExtra
userDefinedResBB = True
xMinResUserDefined = -100
xMaxResUserDefined = 100
yMinResUserDefined = -45
yMaxResUserDefined = 45
zMinResUserDefined = 0
zMaxResUserDefined = zGridSize


nNF = 3      # number of natural fractures
fracHeightRange = [1,2]
angleNFRange = [-math.pi/3,math.pi/3]    # angle in radian

resLengthX = xMaxResUserDefined - xMinResUserDefined
resLengthY = yMaxResUserDefined - yMinResUserDefined
lenNFRange = [10,min(fracLengthScaling*resLengthX,fracLengthScaling*resLengthY)]
mainFractureLength = 40  # main HF length

## mesh writing to openfoam inputs
debug=1      # Print Verbosity (0=silent => 3=chatty)
verify=False # Verify face order, might take longer time

###################################################
###             NF generation                   ###
###################################################

cx = []
cy = []
cz = []
lenNF = []
theta = []
fracHeight = []

for i in range(nNF):
    cx.append(round(random.uniform(-fracCenterExtent*resLengthX/2,fracCenterExtent*resLengthX/2),2))
    cy.append(round(random.uniform(-fracCenterExtent*resLengthY/2,fracCenterExtent*resLengthY/2),2))
    cz.append(round(random.uniform(0,0),2))
    theta.append(round(random.uniform(angleNFRange[0],angleNFRange[1]),2))
    lenNF.append(round(random.uniform(lenNFRange[0],lenNFRange[1]),2))
    fracHeight.append(round(random.uniform(fracHeightRange[0],fracHeightRange[1]),0))

###      If reading NFs from a file            ###
fileNameNF = 'meshingInput.txt'
cx = []
cy = []
cz = []
lenNF = []
theta = []
fracHeight = []
with open(fileNameNF, 'r') as filehandle:
    for line in filehandle:
        if not line.isspace():
            lineText = line.split(",")
            cx.append(float(lineText[0]))
            cy.append(float(lineText[1]))
            cz.append(float(lineText[2]))
            theta.append(float(lineText[3]))
            lenNF.append(float(lineText[4]))
            fracHeight.append(float(lineText[5]))
nNF = len(cx)

# generating co-ordinates to be used in salome geometry
pxA = []
pyA = []
pzA = []
pxB = []
pyB = []
pzB = []

for i in range(nNF):
    pxA.append(cx[i] + 0.5*lenNF[i]*math.cos(theta[i]))
    pyA.append(cy[i] + 0.5*lenNF[i]*math.sin(theta[i]))
    pzA.append(cz[i])
    pxB.append(cx[i] - 0.5*lenNF[i]*math.cos(theta[i]))
    pyB.append(cy[i] - 0.5*lenNF[i]*math.sin(theta[i]))
    pzB.append(cz[i])

###################################################
###             SALOME GEOM                     ###
###################################################

import sys
import salome

salome.salome_init()
print("starting the geometry creation ....")

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder

geompy = geomBuilder.New()

OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

VeticesA = []
VeticesB = []
fracLines = []
fracFaces = []

for j in range(nNF):
    nj = j+1

    PointA = geompy.MakeVertex(pxA[j],pyA[j],pzA[j])
    # geompy.addToStudy(PointA, "VertexA_" + str(nj))
    VeticesA.append(PointA)
    PointB = geompy.MakeVertex(pxB[j],pyB[j],pzB[j])
    # geompy.addToStudy(PointBBottom, "VertexB_" + str(nj))
    VeticesB.append(PointB)

    fracLine = geompy.MakeLineTwoPnt(PointA, PointB)
    # geompy.addToStudy(fracLine, "fracLine_" + str(nj))
    fracLines.append(fracLine)

    fracFace = geompy.MakePrismVecH(fracLine, OZ, fracHeight[j])
    geompy.addToStudy(fracFace, "fracFace_" + str(nj))
    fracFaces.append(fracFace)

fracFaces = geompy.MakeGlueEdges(fracFaces,glueTol)
crackFaces = geompy.MakeCompound(fracFaces)
geompy.addToStudy( crackFaces, 'crackFaces' )
crackFacesForCut = geompy.MakePartition([crackFaces], [], [], [], geompy.ShapeType["FACE"], 0, [], 0)
geompy.addToStudy( crackFacesForCut, 'crackFacesForCut' )

if(userDefinedResBB):
    xMinResBB = xMinResUserDefined
    xMaxResBB = xMaxResUserDefined
    yMinResBB = yMinResUserDefined
    yMaxResBB = yMaxResUserDefined
    zMinResBB = zMinResUserDefined
    zMaxResBB = zMaxResUserDefined

print("xMinResBB = ",xMinResBB)
print("xMaxResBB = ",xMaxResBB)
print("yMinResBB = ",yMinResBB)
print("yMaxResBB = ",yMaxResBB)
print("zMinResBB = ",zMinResBB)
print("zMaxResBB = ",zMaxResBB)

vertex0ResBB = geompy.MakeVertex(xMinResBB,yMinResBB,zMinResBB)
vertex1ResBB = geompy.MakeVertex(xMaxResBB,yMinResBB,zMinResBB)
vertex2ResBB = geompy.MakeVertex(xMaxResBB,yMaxResBB,zMinResBB)
vertex3ResBB = geompy.MakeVertex(xMinResBB,yMaxResBB,zMinResBB)
vertex4ResBB = geompy.MakeVertex(xMinResBB,yMinResBB,zMaxResBB)
vertex5ResBB = geompy.MakeVertex(xMaxResBB,yMinResBB,zMaxResBB)
vertex6ResBB = geompy.MakeVertex(xMaxResBB,yMaxResBB,zMaxResBB)
vertex7ResBB = geompy.MakeVertex(xMinResBB,yMaxResBB,zMaxResBB)

reservoirBox = geompy.MakeBoxTwoPnt(vertex0ResBB, vertex6ResBB)
geompy.addToStudy( reservoirBox, 'reservoirBox' )

xNegEdge = geompy.MakeLineTwoPnt(vertex0ResBB, vertex3ResBB)
xPosEdge = geompy.MakeLineTwoPnt(vertex1ResBB, vertex2ResBB)
yNegEdge = geompy.MakeLineTwoPnt(vertex0ResBB, vertex1ResBB)
yPosEdge = geompy.MakeLineTwoPnt(vertex2ResBB, vertex3ResBB)
topEdge = geompy.MakeLineTwoPnt(vertex4ResBB, vertex5ResBB)
bottomEdge = geompy.MakeLineTwoPnt(vertex0ResBB, vertex1ResBB)

xNegFace = geompy.MakePrismVecH(xNegEdge, OZ, zMaxResBB-zMinResBB)
xPosFace = geompy.MakePrismVecH(xPosEdge, OZ, zMaxResBB-zMinResBB)
yNegFace = geompy.MakePrismVecH(yNegEdge, OZ, zMaxResBB-zMinResBB)
yPosFace = geompy.MakePrismVecH(yPosEdge, OZ, zMaxResBB-zMinResBB)
resTop = geompy.MakePrismVecH(topEdge, OY, yMaxResBB-yMinResBB)
resBottom = geompy.MakePrismVecH(bottomEdge, OY, yMaxResBB-yMinResBB)

# reservoirSolid for 3D meshing
reservoirSolid = geompy.MakePartition([reservoirBox, crackFaces], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
resFracFaceList = geompy.SubShapeAllIDs(reservoirSolid, geompy.ShapeType["FACE"])
# print("IDs of edges of reservoirBase:", resFracFaceList, "(unsorted)")
resFracFaceList = resFracFaceList[6:]
crackFaceGroup = geompy.CreateGroup(reservoirSolid, geompy.ShapeType["FACE"])
geompy.UnionIDs(crackFaceGroup, resFracFaceList)

# reservoirBase for 2D meshing
reservoirBase = geompy.MakeCutList(resBottom, [crackFacesForCut], True)
resFracEdgeList = geompy.SubShapeAllIDs(reservoirBase, geompy.ShapeType["EDGE"])
# print("IDs of edges of reservoirBase:", resFracEdgeList, "(unsorted)")
crackEdgeList = resFracEdgeList[4:]
crackEdgeGroup = geompy.CreateGroup(reservoirBase, geompy.ShapeType["EDGE"])
geompy.UnionIDs(crackEdgeGroup, crackEdgeList)

geompy.addToStudy( crackFaces, 'crackFaces' )
geompy.addToStudy( reservoirSolid, 'reservoirSolid' )
# geompy.addToStudyInFather( reservoirSolid, crackFaces, 'crackFaces' )
geompy.addToStudyInFather( reservoirSolid, crackFaceGroup, 'crackFaceGroup' )
geompy.addToStudyInFather( reservoirSolid, resBottom, 'resBottom' )
geompy.addToStudyInFather( reservoirSolid, resTop, 'resTop' )
geompy.addToStudyInFather( reservoirSolid, xNegFace, 'xNegFace' )
geompy.addToStudyInFather( reservoirSolid, yNegFace, 'yNegFace' )
geompy.addToStudyInFather( reservoirSolid, xPosFace, 'xPosFace' )
geompy.addToStudyInFather( reservoirSolid, yPosFace, 'yPosFace' )

geompy.addToStudy( reservoirBase, 'reservoirBase' )
geompy.addToStudyInFather( reservoirBase, crackEdgeGroup, 'crackEdgeGroup' )
geompy.addToStudyInFather( reservoirBase, resBottom, 'resBottom' )
geompy.addToStudyInFather( reservoirBase, resTop, 'resTop' )
geompy.addToStudyInFather( reservoirBase, xNegFace, 'xNegFace' )
geompy.addToStudyInFather( reservoirBase, yNegFace, 'yNegFace' )
geompy.addToStudyInFather( reservoirBase, xPosFace, 'xPosFace' )
geompy.addToStudyInFather( reservoirBase, yPosFace, 'yPosFace' )

fnameVTK = os.getcwd()+"/crackFaces.vtk"
geompy.ExportVTK(crackFaces, fnameVTK , 0.001)

print('geometry module DONE!!!!!!')


###################################################
###             SALOME MESHING                  ###
###################################################

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

if(isMesh3D):    
    resMesh = smesh.Mesh(reservoirSolid)
    NETGEN_1D_2D_3D = resMesh.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
    NETGEN_3D_Parameters_1.SetMaxSize( resOuterGridSize )
    NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
    NETGEN_3D_Parameters_1.SetOptimize( 1 )
    NETGEN_3D_Parameters_1.SetFineness( SetFineness3D )
    if(SetFineness3D == 5):
        NETGEN_3D_Parameters_1.SetGrowthRate( SetGrowthRate3D )
        NETGEN_3D_Parameters_1.SetNbSegPerEdge( SetNbSegPerEdge3D )
        NETGEN_3D_Parameters_1.SetNbSegPerRadius( SetNbSegPerRadius3D )
    NETGEN_3D_Parameters_1.SetChordalError( 0.1 )
    NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
    NETGEN_3D_Parameters_1.SetMinSize( resGridMinSize )
    NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
    NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
    NETGEN_3D_Parameters_1.SetQuadAllowed( useQuadElementsRes )
    Regular_1D = resMesh.Segment(geom=crackFaceGroup)
    Local_Length_1 = Regular_1D.LocalLength(fracGridMinSize,None,1e-07)
    NETGEN_2D = resMesh.Triangle(algo=smeshBuilder.NETGEN_2D,geom=crackFaceGroup)
    NETGEN_2D_Parameters_1 = NETGEN_2D.Parameters()
    NETGEN_2D_Parameters_1.SetMaxSize( fracGridMaxSize )
    NETGEN_2D_Parameters_1.SetOptimize( 1 )
    NETGEN_2D_Parameters_1.SetFineness( 2 )
    NETGEN_2D_Parameters_1.SetChordalError( 0.1 )
    NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
    NETGEN_2D_Parameters_1.SetMinSize( fracGridMinSize )
    NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
    NETGEN_2D_Parameters_1.SetQuadAllowed( useQuadElementsFrac )
    isDone = resMesh.Compute()
    crackFaceMesh = Regular_1D.GetSubMesh()
else:
    resMesh = smesh.Mesh(reservoirBase)
    Regular_1D = resMesh.Segment()
    Local_Length_2 = Regular_1D.LocalLength(resOuterGridSize,None,1e-007)

    NETGEN_2D = resMesh.Triangle(algo=smeshBuilder.NETGEN_2D)

    NETGEN_2D_Parameters_1 = NETGEN_2D.Parameters()
    NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
    NETGEN_2D_Parameters_1.SetOptimize( 1 )
    NETGEN_2D_Parameters_1.SetFineness( SetFineness2D )
    if(SetFineness2D == 5):
        NETGEN_2D_Parameters_1.SetGrowthRate( SetGrowthRate2D )
        NETGEN_2D_Parameters_1.SetNbSegPerEdge( SetNbSegPerEdge2D )
        NETGEN_2D_Parameters_1.SetNbSegPerRadius( SetNbSegPerRadius2D )
    NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
    NETGEN_2D_Parameters_1.SetQuadAllowed( useQuadElementsFrac )

    Regular_1D = resMesh.Segment(geom=crackEdgeGroup)
    Local_Length_1 = Regular_1D.LocalLength(fracGridMinSize,None,1e-007)
    NETGEN_2D_Parameters_1.SetMaxSize( maxGridSize )
    NETGEN_2D_Parameters_1.SetMinSize( minGridSize )
    NETGEN_2D_Parameters_1.SetFuseEdges( 184 )
    isDone = resMesh.Compute()
    resMesh.ExtrusionSweepObjects( [ resMesh ], [ resMesh ], [ resMesh ], [ 0, 0, zGridSize ], timesExtrusion, 1 )

    subMeshCrack = Regular_1D.GetSubMesh()

# defining crack faces
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,crackFaces,SMESH.FT_Undefined,SMESH.FT_LogicalAND,1e-002)
aCriteria.append(aCriterion)
# aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_EntityType,SMESH.FT_Undefined,SMESH.Entity_Quadrangle,SMESH.FT_Undefined,SMESH.FT_LogicalAND,1e-007)
# aCriteria.append(aCriterion)
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,resTop,SMESH.FT_LogicalNOT,SMESH.FT_LogicalAND,1e-009)
aCriteria.append(aCriterion)
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,resBottom,SMESH.FT_LogicalNOT,SMESH.FT_LogicalAND,1e-009)
aCriteria.append(aCriterion)
aFilterCrack = smesh.GetFilterFromCriteria(aCriteria)
aFilterCrack.SetMesh(resMesh.GetMesh())
crack = resMesh.GroupOnFilter( SMESH.FACE, 'crack', aFilterCrack )
crack.SetName( 'crack' )

# defining topBottom
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,resBottom,SMESH.FT_Undefined,SMESH.FT_LogicalOR,1e-009)
aCriteria.append(aCriterion)
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,resTop,SMESH.FT_Undefined,SMESH.FT_Undefined,1e-009)
aCriteria.append(aCriterion)
aFilterResTopBottom = smesh.GetFilterFromCriteria(aCriteria)
aFilterResTopBottom.SetMesh(resMesh.GetMesh())
topBottom = resMesh.GroupOnFilter( SMESH.FACE, 'topBottom', aFilterResTopBottom )
topBottom.SetName( 'topBottom' )

# defining xNegPos
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,xNegFace,SMESH.FT_Undefined,SMESH.FT_LogicalOR,1e-009)
aCriteria.append(aCriterion)
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,xPosFace,SMESH.FT_Undefined,SMESH.FT_Undefined,1e-009)
aCriteria.append(aCriterion)
aFilterXNegPos = smesh.GetFilterFromCriteria(aCriteria)
aFilterXNegPos.SetMesh(resMesh.GetMesh())
xNegPos = resMesh.GroupOnFilter( SMESH.FACE, 'xNegPos', aFilterXNegPos )
xNegPos.SetName( 'xNegPos' )


# defining yNegPos
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,yNegFace,SMESH.FT_Undefined,SMESH.FT_LogicalOR,1e-009)
aCriteria.append(aCriterion)
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,yPosFace,SMESH.FT_Undefined,SMESH.FT_Undefined,1e-009)
aCriteria.append(aCriterion)
aFilterYNegPos = smesh.GetFilterFromCriteria(aCriteria)
aFilterYNegPos.SetMesh(resMesh.GetMesh())
yNegPos = resMesh.GroupOnFilter( SMESH.FACE, 'yNegPos', aFilterYNegPos )
yNegPos.SetName( 'yNegPos' )

if(isMesh3D):    
    ## Set names of Mesh objects
    smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
    smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
    smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
    smesh.SetName(Local_Length_1, 'Local Length_1')
    smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
    smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
    smesh.SetName(resMesh.GetMesh(), 'resMesh')
    smesh.SetName(crackFaceMesh, 'crackFaceMesh')
else:
    ## Set names of Mesh objects
    smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
    smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
    smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
    smesh.SetName(Local_Length_2, 'Local Length_2')
    smesh.SetName(Local_Length_1, 'Local Length_1')
    smesh.SetName(resMesh.GetMesh(), 'resMesh')
    smesh.SetName(crack, 'crack')
    smesh.SetName(subMeshCrack, 'subMeshCrack')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()


# ###################################################
# ##  Functions to WRITE Mesh in OpenFOAM Format   ##
# ###################################################

import sys
import salome
import SMESH
from salome.smesh import smeshBuilder
import os, time

class MeshBuffer(object):
    """
    Limits the calls to Salome by buffering the face and key details of volumes to speed up exporting
    """

    def __init__(self,mesh,v):
        i=0
        faces=list()
        keys=list()
        fnodes=mesh.GetElemFaceNodes(v,i)

        while fnodes:                           #While not empty list
            faces.append(fnodes)                #Face list
            keys.append(tuple(sorted(fnodes)))  #Buffer key
            i += 1
            fnodes=mesh.GetElemFaceNodes(v,i)

        self.v=v         #The volume
        self.faces=faces #The sorted face list
        self.keys=keys
        self.fL=i        #The number of faces

    @staticmethod
    def Key(fnodes):
        """Takes the nodes and compresses them into a hashable key"""
        return tuple(sorted(fnodes))


    @staticmethod
    def ReverseKey(fnodes):
        """Takes the nodes and compresses them into a hashable key reversed for baffles"""
        if(type(fnodes) is tuple): return tuple(reversed(fnodes))
        else: return tuple(sorted(fnodes, reverse=True))

def exportToFoam(mesh,dirname='polyMesh'):
    """
    Export a mesh to OpenFOAM.

    args:
        +    mesh: The mesh
        + dirname: The mesh directory to write to

    The algorithm works as follows:
    [1] Loop through the boundaries and collect all faces in each group.
        Faces that don't have a group will be added to the group defaultPatches.

    [2] Loop through all cells (volumes) and each face in the cell.
        If the face has been visited before it is added to the neighbour list
        If not, then check it is a boundary face.
            If it is, add the cell to the end of owner.
            If not a boundary face and has not yet been visited add it to the list of internal faces.

    To check if a face has been visited a dictionary is used.
    The key is the sorted list of face nodes converted to a string.
    The value is the face id. Eg: facesSorted[key] = value
    """

    starttime=time.time()

    #try to open files
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    try:
        filePoints=open(dirname +"/_points_",'w')
        fileFaces=open(dirname +"/_faces_",'w')
        fileOwner=open(dirname + "/_owner_",'w')
        fileNeighbour=open(dirname + "/_neighbour_",'w')
        fileBoundary=open(dirname + "/_boundary_",'w')
        fileCrackInfo=open(os.getcwd() + "/crackInfo",'w')
    except Exception:
        print("could not open files aborting")
        return

    #Get salome properties
    smesh = smeshBuilder.New()
    __debugPrint__('Number of nodes: %d\n' %(mesh.NbNodes()))

    volumes=mesh.GetElementsByType(SMESH.VOLUME)
    __debugPrint__("Number of cells: %d\n" %len(volumes))

    __debugPrint__('Counting number of faces:\n')
    #Filter faces
    filter=smesh.GetFilter(SMESH.EDGE,SMESH.FT_FreeFaces)
    extFaces=mesh.GetIdsFromFilter(filter)
    nrBCfaces=len(extFaces)
    nrExtFaces=len(extFaces)
    #nrBCfaces=mesh.NbFaces();#number of bcfaces in Salome

    nrFaces=0
    buffers=list()
    for v in volumes:
        b = MeshBuffer(mesh, v)
        nrFaces += b.fL
        buffers.append(b)

     #all internal faces will be counted twice, external faces once
    #so:
    nrFaces=(nrFaces+nrExtFaces)/2
    nrIntFaces=nrFaces-nrBCfaces #
    __debugPrint__('total number of faces: %d, internal: %d, external %d\n'  \
        %(nrFaces,nrIntFaces,nrExtFaces))

    __debugPrint__("Converting mesh to OpenFOAM\n")
    faces=[] #list of internal face nodes ((1 2 3 4 ... ))
    facesSorted=dict() #each list of nodes is sorted.
    bcFaces=[] #list of bc faces (
    bcFacesSorted=dict()
    owner=[] #owner file, (of face id, volume id)
    neighbour=[] #neighbour file (of face id, volume id) only internal faces

    #Loop over all salome boundary elemets (faces)
    # and store them inte the list bcFaces
    grpStartFace=[] # list of face ids where the BCs starts
    grpNrFaces=[] # list of number faces in each BC
    grpNames=[] #list of the group name.
    ofbcfid=0   # bc face id in openfoam
    nrExtFacesInGroups=0

    for gr in mesh.GetGroups():
        if gr.GetType() == SMESH.FACE:
            grpNames.append(gr.GetName())

            # mesh for fracture propagation
            if ("crack" in gr.GetName() and dynamicCrack):
                # print('found group \"%s\" of type %s, %d\n' \
                #                %(gr.GetName(),gr.GetType(),len(gr.GetIDs())))
                grIds=gr.GetIDs()
                #loop over faces in group
                for sfid in grIds:
                    fnodes=mesh.GetElemNodes(sfid)
                    # write crackInfo if group name is crack
                    for crackPI in fnodes:
                        pos=mesh.GetNodeXYZ(crackPI)
                        fileCrackInfo.write("%g, %g, %g" %(pos[0],pos[1],pos[2]))
                        if(crackPI != fnodes[-1]):
                            fileCrackInfo.write(";")
                    fileCrackInfo.write("\n")
            # static fractures
            else:
                # __debugPrint__('found group \"%s\" of type %s, %d\n' \
                #                %(gr.GetName(),gr.GetType(),len(gr.GetIDs())),2)
                print('found group \"%s\" of type %s, %d\n' \
                               %(gr.GetName(),gr.GetType(),len(gr.GetIDs())))
                grIds=gr.GetIDs()
                nr=len(grIds)
                if  nr > 0:
                    grpStartFace.append(nrIntFaces+ofbcfid)
                    grpNrFaces.append(nr)

                #loop over faces in group

                for sfid in grIds:
                    fnodes=mesh.GetElemNodes(sfid)
                    key=MeshBuffer.Key(fnodes)
                    if not key in bcFacesSorted:
                        bcFaces.append(fnodes)
                        bcFacesSorted[key]=ofbcfid
                        ofbcfid+=1
                    else:
                        raise Exception(\
                            "Error the face, elemId %d, %s belongs to two " %(sfid,fnodes)  +\
                                "or more groups. One is : %s"  %(gr.GetName()))

                    # write crackInfo if group name is crack
                    if ("crack" in gr.GetName()):
                        # fileCrackInfo.write("[")
                        for crackPI in fnodes:
                            pos=mesh.GetNodeXYZ(crackPI)
                            fileCrackInfo.write("%g, %g, %g" %(pos[0],pos[1],pos[2]))
                            if(crackPI != fnodes[-1]):
                                fileCrackInfo.write(";")
                        fileCrackInfo.write("\n")

                #if the group is a baffle then the faces should be added twice
                if __isGroupBaffle__(mesh,gr,extFaces):
                    nrBCfaces+=nr
                    nrFaces+=nr
                    nrIntFaces-=nr
                    #since nrIntFaces is reduced all previously grpStartFaces are
                    #out of sync
                    grpStartFace=[x-nr for x in grpStartFace]
                    grpNrFaces[-1]=nr*2
                    for sfid in gr.GetIDs():
                        fnodes=mesh.GetElemNodes(sfid)
                        key=MeshBuffer.ReverseKey(fnodes)
                        bcFaces.append(fnodes)
                        bcFacesSorted[key]=ofbcfid
                        ofbcfid+=1
                else:
                    nrExtFacesInGroups+=nr

    __debugPrint__('total number of faces: %d, internal: %d, external %d\n'  \
        %(nrFaces,nrIntFaces,nrExtFaces),2)

    #Do the defined groups cover all BC-faces?
    if nrExtFacesInGroups < nrExtFaces:
        __debugPrint__("Warning, some elements don't have a group (BC). " +\
                       "Adding to a new group called defaultPatches\n",1)

        grpStartFace.append(nrIntFaces+ofbcfid)
        grpNrFaces.append(nrExtFaces-nrExtFacesInGroups)

        salomeIDs=[]
        for face in extFaces:
            fnodes=mesh.GetElemNodes(face)
            key=MeshBuffer.Key(fnodes)
            try:
                bcFacesSorted[key]
            except KeyError:
                #if not in dict then add to default patches
                bcFaces.append(fnodes)
                bcFacesSorted[key]=ofbcfid
                salomeIDs.append(face)
                ofbcfid+=1

        newGrpName="defaultPatches"
        nri=1
        while newGrpName in grpNames:
            newGrpName="defaultPatches_%d" %nri
            nri+=1

        grpNames.append(newGrpName)
        #function might have different name
        try:
            defGroup=mesh.CreateGroup(SMESH.FACE, newGrpName )
        except AttributeError:
            defGroup=mesh.CreateEmptyGroup(SMESH.FACE, newGrpName )

        defGroup.Add(salomeIDs)
        smesh.SetName(defGroup, newGrpName)
        if salome.sg.hasDesktop():
            salome.sg.updateObjBrowser()

    #initialise the list faces vs owner/neighbour cells
    nrFaces = round(nrFaces)
    nrIntFaces = round(nrIntFaces)
    owner=[-1]*nrFaces
    neighbour=[-1]*nrIntFaces

    __debugPrint__("Finished processing boundary faces\n")
    __debugPrint__('bcFaces: %d\n' %(len(bcFaces)),2)
    __debugPrint__(str(bcFaces)+"\n",3)
    __debugPrint__('bcFacesSorted: %d\n' %(len(bcFacesSorted)),2)
    __debugPrint__(str(bcFacesSorted)+"\n",3)
    __debugPrint__('owner: %d\n' %(len(owner)),2)
    __debugPrint__(str(owner)+"\n",3)
    __debugPrint__('neighbour: %d\n' %(len(neighbour)),2)
    __debugPrint__(str(neighbour)+"\n",3)

    offid=0;
    ofvid=0; #volume id in openfoam
    for b in buffers:
        if debug > 2: #Salome call only if verbose
            nodes=mesh.GetElemNodes(b.v)
            __debugPrint__('volume id: %d, num nodes %d, nodes:%s \n' %(b.v,len(nodes),nodes),3)

        fi = 0 #Face index
        while fi < b.fL:
            fnodes=b.faces[fi]
            key=b.keys[fi]
            #Check if the node is already in list
            try:
                fidinof=facesSorted[key]
                #if faceSorted didn't throw an exception then the face is
                #already in the dict. Its an internal face and should be added
                # to the neighbour list
                #print "fidinof %d" %fidinof

                neighbour[fidinof]=ofvid
                __debugPrint__('\tan owner already exist for %d, %s, cell %d\n' %(fi,fnodes,ofvid),3)
            except KeyError:
                #the face is not in the list of internal faces
                #it might a new face or a BCface.
                try:
                    bcind=bcFacesSorted[key]
                    #if no exception was trown then it's a bc face
                    __debugPrint__('\t found bc face: %d, %s, cell %d\n' %(bcind,fnodes,ofvid),3)
                    #if the face belongs to a baffle then it exits twice in owner
                    #check dont overwrite owner
                    if owner[nrIntFaces+bcind]==-1:
                        owner[nrIntFaces+bcind]=ofvid
                        bcFaces[bcind]=fnodes
                    else:
                        #build functions that looks for baffles in bclist. with bcind
                        key=MeshBuffer.ReverseKey(fnodes)
                        bcind=bcFacesSorted[key]
                        #make sure the faces has the correct orientation
                        bcFaces[bcind]=fnodes
                        owner[nrIntFaces+bcind]=ofvid

                except KeyError:
                    #the face is not in bc list either so it's a new internal face
                    __debugPrint__('\t a new face was found, %d, %s, cell %d\n' %(fi,fnodes,ofvid),3)
                    if verify:
                        if not __verifyFaceOrder__(mesh,nodes,fnodes):
                            __debugPrint__("\t face has bad order, reversing order\n",3)
                            fnodes.reverse()

                    faces.append(fnodes)
                    key=b.keys[fi]
                    facesSorted[key]=offid
                    owner[offid]=ofvid

                    offid+=1
                    if(nrFaces > 50 and offid % (nrFaces/50)==0):
                        if(offid % ((nrFaces/50)*10) == 0):
                            __debugPrint__(':',1)
                        else:
                            __debugPrint__('.',1)
            fi+=1
        ofvid+=1

        # end for v in volumes

    nrCells=ofvid
    __debugPrint__("Finished processing volumes.\n")
    __debugPrint__('faces: %d\n' %(len(faces)),2)
    __debugPrint__(str(faces)+"\n",3)
    __debugPrint__('facesSorted: %d\n' %(len(facesSorted)),2)
    __debugPrint__(str(facesSorted)+"\n",3)
    __debugPrint__('owner: %d\n' %(len(owner)),2)
    __debugPrint__(str(owner)+"\n",3)
    __debugPrint__('neighbour: %d\n' %(len(neighbour)),2)
    __debugPrint__(str(neighbour)+"\n",3)

    #Convert to "upper triangular order"
    #owner is sorted, for each cell sort faces it's neighbour faces
    # i.e. change
    # owner   neighbour          owner   neighbour
    #     0          15                    0                3
    #     0            3          to       0              15
    #     0           17                   0              17
    #     1            5                    1                5
    # any changes made to neighbour are repeated to faces.

    __debugPrint__("Sorting faces in upper triangular order\n",1)
    ownedfaces=1
    for faceId in range(0,nrIntFaces):
        cellId=owner[faceId]
        nextCellId=owner[faceId+1] #np since len(owner) > nrIntFaces
        if cellId == nextCellId:
            ownedfaces+=1
            continue

        if ownedfaces > 1:
            sId=faceId-ownedfaces+1 #start ID
            eId=faceId #end ID
            inds=list(range(sId,eId+1))
            inds.sort(key=neighbour.__getitem__)
            neighbour[sId:eId+1]=map(neighbour.__getitem__,inds)
            faces[sId:eId+1]=map(faces.__getitem__,inds)

        ownedfaces=1
    converttime=time.time()-starttime

    #WRITE points to file
    __debugPrint__("Writing the file points\n")
    __writeHeader__(filePoints,"points")
    points=mesh.GetElementsByType(SMESH.NODE)
    nrPoints=len(points)
    filePoints.write("\n%d\n(\n" %(nrPoints))

    for n,ni in enumerate(points):
        pos=mesh.GetNodeXYZ(ni)
        filePoints.write("\t(%g %g %g)\n" %(pos[0],pos[1],pos[2]))

    filePoints.write(")\n")
    filePoints.flush()
    filePoints.close()

    #WRITE faces to file
    __debugPrint__("Writing the file faces\n")
    __writeHeader__(fileFaces,"faces")
    fileFaces.write("\n%d\n(\n" %(nrFaces))

    for node in faces:
        fileFaces.write("\t%d(" %(len(node)))
        for p in node:
            #salome starts to count from one, OpenFOAM from zero
            fileFaces.write("%d " %(p-1))
        fileFaces.write(")\n")

    #internal nodes are done output bcnodes
    for node in bcFaces:
        fileFaces.write("\t%d(" %(len(node)))
        for p in node:
            #salome starts to count from one, OpenFOAM from zero
            fileFaces.write("%d " %(p-1))
        fileFaces.write(")\n")

    fileFaces.write(")\n")
    fileFaces.flush()
    fileFaces.close()

    #WRITE owner to file
    __debugPrint__("Writing the file owner\n")
    __writeHeader__(fileOwner,"owner",nrPoints,nrCells,nrFaces,nrIntFaces)
    fileOwner.write("\n%d\n(\n" %(len(owner)))

    for cell in owner:
        fileOwner.write(" %d \n" %(cell))
    fileOwner.write(")\n")
    fileOwner.flush()
    fileOwner.close()

    #WRITE neighbour
    __debugPrint__("Writing the file neighbour\n")
    __writeHeader__(fileNeighbour,"neighbour",nrPoints,nrCells,nrFaces,nrIntFaces)

    fileNeighbour.write("\n%d\n(\n" %(len(neighbour)))
    for cell in neighbour:
        fileNeighbour.write(" %d\n" %(cell))

    fileNeighbour.write(")\n")
    fileNeighbour.flush()
    fileNeighbour.close()

    #WRITE boundary file
    __debugPrint__("Writing the file boundary\n")
    __writeHeader__(fileBoundary,"boundary")

    fileBoundary.write("%d\n(\n" %len(grpStartFace))
    ind = -1
    for gname in grpNames:
        # mesh for fracture propagation
        if ("crack" in gname and dynamicCrack):
            # do nothing
            print("\n\ncrack not created in the boundary file\n\n")
        else:
            ind = ind + 1
            fileBoundary.write("\t%s\n\t{\n" %gname)
            fileBoundary.write("\ttype\t\t")
            if "wall" in gname.lower():
                fileBoundary.write("wall;\n")
            elif "crack" in gname.lower():
                fileBoundary.write("cohesiveRegionCouple;\n")
            elif "topbottom" in gname.lower():
                if timesExtrusion == 1:
                    fileBoundary.write("empty;\n")
                else:
                    fileBoundary.write("patch;\n")
            else:
                fileBoundary.write("patch;\n")
            fileBoundary.write("\tnFaces\t\t%d;\n" %grpNrFaces[ind])
            fileBoundary.write("\tstartFace\t%d;\n" %grpStartFace[ind])
            fileBoundary.write("\t}\n")
    fileBoundary.write(")\n")
    fileBoundary.close()

    #WRITE cellZones
#Count number of cellZones
    nrCellZones=0;
    cellZonesName=list();
    for grp in mesh.GetGroups():
        if grp.GetType() == SMESH.VOLUME:
            nrCellZones+=1
            cellZonesName.append(grp.GetName())

    if nrCellZones > 0:
        try:
            fileCellZones=open(dirname + "/cellZones",'w')
        except Exception:
            print("Could not open the file cellZones, other files are ok.")
        __debugPrint__("Writing file cellZones\n")
        #create a dictionary where salomeIDs are keys
        #and OF cell ids are values.
        scToOFc=dict([sa,of] for of,sa in enumerate(volumes))
        __writeHeader__(fileCellZones,"cellZones")
        fileCellZones.write("\n%d(\n" %nrCellZones)

        for grp in mesh.GetGroups():
            if grp.GetType() == SMESH.VOLUME:
                fileCellZones.write(grp.GetName()+"\n{\n")
                fileCellZones.write("\ttype\tcellZone;\n")
                fileCellZones.write("\tcellLabels\tList<label>\n")
                cellSalomeIDs=grp.GetIDs()
                nrGrpCells=len(cellSalomeIDs)
                fileCellZones.write("%d\n(\n" %nrGrpCells)

                for csId in cellSalomeIDs:
                    ofID=scToOFc[csId]
                    fileCellZones.write("%d\n" %ofID)

                fileCellZones.write(");\n}\n")

        fileCellZones.write(")\n")
        fileCellZones.flush()
        fileCellZones.close()

    totaltime=time.time()-starttime
    __debugPrint__("Finished writing to %s/%s \n" %(os.getcwd(),dirname))
    __debugPrint__("Converted mesh in %.0fs\n" %(converttime),1)
    __debugPrint__("Wrote mesh in %.0fs\n" %(totaltime-converttime),1)
    __debugPrint__("Total time: %0.fs\n" %totaltime,1)


def __writeHeader__(file,fileType,nrPoints=0,nrCells=0,nrFaces=0,nrIntFaces=0):
    """Write a header for the files points, faces, owner, neighbour"""

    file.write("/*" + "-"*68 + "*\\\n" )
    file.write("|" + " "*70 + "|\n")
    file.write("|" + " "*4 + "File exported from Salome Platform" +\
                   " using SalomeToFoamExporter" +" "*5 +"|\n")
    file.write("|" + " "*70 + "|\n")
    file.write("\*" + "-"*68 + "*/\n")

    file.write("FoamFile\n{\n")
    file.write("\tversion\t\t2.0;\n")
    file.write("\tformat\t\tascii;\n")
    file.write("\tclass\t\t")

    if(fileType =="points"):
        file.write("vectorField;\n")
    elif(fileType =="faces"):
        file.write("faceList;\n")
    elif(fileType =="owner" or fileType=="neighbour"):
        file.write("labelList;\n")
        file.write("\tnote\t\t\"nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d\";\n" \
                       %(nrPoints,nrCells,nrFaces,nrIntFaces))

    elif(fileType == "boundary"):
        file.write("polyBoundaryMesh;\n")
    elif(fileType=="cellZones"):
        file.write("regIOobject;\n")

    file.write("\tlocation\t\"constant/polyMesh\";\n")
    file.write("\tobject\t\t" + fileType +";\n")
    file.write("}\n\n")

def __debugPrint__(msg,level=1):
    """Print only if level >= debug """
    if(debug >= level ):
        print(msg,)

def __verifyFaceOrder__(mesh,vnodes,fnodes):
    """
    Verify if the face order is correct. I.e. pointing out of the cell
    calc vol center
    calc f center
    calc ftov=fcenter-vcenter
     calc fnormal=first to second cross first to last
    if ftov dot fnormal >0 reverse order

    """
    vc=__cog__(mesh,vnodes)
    fc=__cog__(mesh,fnodes)
    fcTovc=__diff__(vc,fc)
    fn=__calcNormal__(mesh,fnodes)
    if(__dotprod__(fn,fcTovc)>0.0):
        return False
    else:
        return True


def __cog__(mesh,nodes):
    """
    calculate the center of gravity.
    """
    c=[0.0,0.0,0.0]
    for n in nodes:
        pos=mesh.GetNodeXYZ(n)
        c[0]+=pos[0]
        c[1]+=pos[1]
        c[2]+=pos[2]

    c[0]/=len(nodes)
    c[1]/=len(nodes)
    c[2]/=len(nodes)
    return c

def __calcNormal__(mesh,nodes):
    """
    Calculate and return face normal.
    """
    p0=mesh.GetNodeXYZ(nodes[0])
    p1=mesh.GetNodeXYZ(nodes[1])
    pn=mesh.GetNodeXYZ(nodes[-1])
    u=__diff__(p1,p0)
    v=__diff__(pn,p0)
    return __crossprod__(u,v)



def __diff__(u,v):
    """
    u - v, in 3D
    """
    res=[0.0]*3
    res[0]=u[0]-v[0]
    res[1]=u[1]-v[1]
    res[2]=u[2]-v[2]
    return res

def __dotprod__(u,v):
    """
    3D scalar dot product
    """
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def __crossprod__(u,v):
    """
    3D cross product
    """
    res=[0.0]*3
    res[0]=u[1]*v[2]-u[2]*v[1]
    res[1]=u[2]*v[0]-u[0]*v[2]
    res[2]=u[0]*v[1]-u[1]*v[0]
    return res

def findSelectedMeshes():
    meshes=list()
    smesh = smeshBuilder.New()
    nrSelected=salome.sg.SelectedCount() # total number of selected items

    foundMesh=False
    for i in range(nrSelected):
        selected=salome.sg.getSelected(i)
        selobjID=salome.myStudy.FindObjectID(selected)
        selobj=selobjID.GetObject()

        if selobj.__class__ == SMESH._objref_SMESH_Mesh or selobj.__class__ == salome.smesh.smeshBuilder.meshProxy :
            mName=selobjID.GetName().replace(" ","_")
            foundMesh=True
            mesh=smesh.Mesh(selobj)
            meshes.append(mesh)

    if not foundMesh:
        print("You have to select a mesh object and then run this script.")
        print("or run the export function directly from TUI")
        print(" import SalomeToOpenFOAM")
        print(" SalomeToOpenFOAM.exportToFoam(mesh,path)")
        return list()
    else:
        return meshes

def __isGroupBaffle__(mesh,group,extFaces):
    for sid in group.GetIDs():
        if not sid in extFaces:
            __debugPrint__("group %s is a baffle\n" %group.GetName(),1)
            return True
    return False


###################################################
###      WRITE Mesh in OpenFOAM Format          ###
###################################################
outdir=os.getcwd()+"/constant/org_polyMesh"
print("outdir of mesh is " + outdir)
# __debugPrint__("found selected mesh exporting to " + outdir + ".\n",1)
exportToFoam(resMesh,outdir)
__debugPrint__("finished exporting\n",2)
#Note: to skip renumberMesh just sort owner
#while moving positions also move neighbour,faces, and bcfaces
#will probably have to first sort the internal faces then bc-faces within each bc

import os
from killSalomeWithPort import killMyPort
killMyPort(os.getenv('NSPORT'))