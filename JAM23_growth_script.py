from abaqus import *
from abaqusConstants import *
import __main__
import part
import material
import section
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import connectorBehavior
import time
import string
import random
import math
import numpy as np
import xyPlot
import displayGroupOdbToolset as dgo
import regionToolset
import displayGroupMdbToolset as dgm
import connectorBehavior

# Units: m, Pa, N,

#######################################################################################
#######################################################################################
def Create_Bilayered_Rectangle(ModelName, PartName, Dimensions):
    cliCommand("""session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)""")

    # Define Model
    mdb.Model(name=ModelName, modelType=STANDARD_EXPLICIT)    
    #mdb.models.changeKey(fromName='Model-1', toName=ModelName)
    s = mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=100.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)

    # Sketch Rectangle
    s.rectangle(point1=(0.0, 0.0), point2=(Length, -Width)) 
    p = mdb.models[ModelName].Part(name=PartName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=Depth)
    
    # Partition Rectangle
    e = p.edges
    c = p.cells    
    
    e_1 = e.findAt(((0., -Width/2., 0.),), ((Length, -Width/2., Depth ),))
    p.PartitionEdgeByParam(edges=e_1, parameter=partition)
    
    e_2 = e.findAt(((0., -Width/2., Depth ),), ((Length, -Width/2., 0.),))
    p.PartitionEdgeByParam(edges=e_2, parameter=1-partition)
         
    pickedCells = c.findAt(((0.0, -Width/2., Depth/2), ))
    v, e, d = p.vertices, p.edges, p.datums
    v_br = v.findAt(coordinates=(Length,  -h, 0.0))
    v_fr = v.findAt(coordinates=(Length,  -h, Depth  ))
    v_fl = v.findAt(coordinates=(0.0, -h, Depth  ))
    p.PartitionCellByPlaneThreePoints(point1=v_br, point2=v_fr, point3=v_fl, cells=pickedCells)

#######################################################################################
#######################################################################################
def Create_Rigid_wall(ModelName):
    # Rigid plane
    s1 = mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=100.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)

    s1.rectangle(point1=(0.0, -20.0), point2=(-2.0, 10.0))
    p = mdb.models[ModelName].Part(name='wall', dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
    p = mdb.models[ModelName].parts['wall']
    p.BaseShell(sketch=s1)
    # Reference point
    v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=(-1.0, -10.0, 0.0))

#######################################################################################
#######################################################################################
def Create_Step(ModelName, Step):
    #mdb.models[ModelName].ExplicitDynamicsStep(name=StepName, previous='Initial', timePeriod= TotalTime, improvedDtMethod = ON)
    mdb.models[ModelName].ExplicitDynamicsStep(name=StepName, previous='Initial', timePeriod= TotalTime, massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, MS, 0.0, None, 0, 0, 0.0, 0.0, 0, None),),)
    mdb.models[ModelName].steps[StepName].setValues(timePeriod=TotalTime, scaleFactor=1.0, linearBulkViscosity=0.0, quadBulkViscosity=0.0, improvedDtMethod=ON)
    mdb.models[ModelName].TabularAmplitude(name='Amp-1', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))

#######################################################################################
#######################################################################################
def Create_Material(ModelName, Materials, ratio): 
    mdb.models[ModelName].Material(name='WHIT', description='*****************************************************************\n  Specification Of Material Properties\n*****************************************************************\n\n\nCOMMENTS FROM *USER MATERIAL\n============================\n\n*      mu        = props(1)\n*      lambda    = props(2)\n*      Gaxn      = props(3)\n*      lam0      = props(4)\n*      **gray matter\n*      mu_g      = props(5)\n*      lambda_g  = props(6)\n*      Gctx      = props(7)\n*      Aconst    = props(8)\n*      charlength = props(9)')
    mdb.models[ModelName].materials['WHIT'].Density(table=((density, ), ))
    mdb.models[ModelName].materials['WHIT'].Depvar(n=4)
    mdb.models[ModelName].materials['WHIT'].UserMaterial(mechanicalConstants=(mu_subcortex, lambda_subcortex, 0.0, 1.0, ratio*mu_subcortex, ratio*lambda_subcortex, Gctx, 0.0, 1.5e-3))

    mdb.models[ModelName].Material(name='GRAY', description='COMMENTS FROM *USER MATERIAL*\n============================\n\n*   mu      = prop(1)\n*    lamda   = prop(2)\n*    Gaxn        = prop(3)\n*    lam0        = prop(4)\n*    **gray matter\n*    mu_g        = prop(5)\n*    lamda_g = prop(6)\n*    Gcxn        = prop(7)\n*    Aconst  = prop(8)\n*    charlength = prop(9)    \n')
    mdb.models[ModelName].materials['GRAY'].Density(table=((density, ), ))
    mdb.models[ModelName].materials['GRAY'].Depvar(n=9)
    mdb.models[ModelName].materials['GRAY'].UserMaterial(mechanicalConstants=(mu_subcortex, lambda_subcortex, 0.0, 1.0, ratio*mu_subcortex, ratio*lambda_subcortex, Gctx, 0.0, 1.5e-3))

#######################################################################################
#######################################################################################
def Create_Section(ModelName, PartName, Dimensions, Materials):
    p = mdb.models[ModelName].parts[PartName]
    c = p.cells 

    c_subcort = c.findAt(((0.0, -Width/2., Depth), )) 
    c_cortex  = c.findAt(((0.0, -h/2.,  Depth), ))

    mdb.models[ModelName].HomogeneousSolidSection(name='Subcortex', material=SubCortMaterial, thickness=None)
    region = p.Set(cells=c_subcort, name='Subcortex')
    p.SectionAssignment(region=region, sectionName='Subcortex', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    
    mdb.models[ModelName].HomogeneousSolidSection(name='Cortex', material=CorticalMaterial, thickness=None)
    region = p.Set(cells=c_cortex, name='Cortex')
    p.SectionAssignment(region=region, sectionName='Cortex', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

#######################################################################################
#######################################################################################           
def Create_Assembly(ModelName, PartName, InstanceName, Dimensions):
    p = mdb.models[ModelName].parts[PartName]
    p1 = mdb.models[ModelName].parts['wall']
    a = mdb.models[ModelName].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)

    a.Instance(name=InstanceName, part=p, dependent=ON)
    a.Instance(name='wall-1', part=p1, dependent=ON)
    a.Instance(name='wall-2', part=p1, dependent=ON)

    a.rotate(instanceList=('wall-1', ), axisPoint=(0.0, 0.0, 0.1), axisDirection=(0.0, -2.0, 0.0), angle=-90.0)
    a.rotate(instanceList=('wall-2', ), axisPoint=(0.0, 0.0, 0.1), axisDirection=(0.0, -2.0, 0.0), angle=-90.0)

    a.translate(instanceList=('wall-1', ), vector=(Length+0.1, 0.0, -0.975))
    a.translate(instanceList=('wall-2', ), vector=(0.1, 0.0, -0.975))

####################################################################################### 
####################################################################################### 
def Create_Sets(ModelName, PartName, Dimensions):
    a = mdb.models[ModelName].rootAssembly
    p = mdb.models[ModelName].parts[PartName]
    f = p.faces
    e = p.edges
    c = p.cells
    v = p.vertices

    # Define coordinate of set-node
    verts = v.findAt(((Length, 0.0, Depth), ))

    # Define coordinate of the sets-faces 
    f_back = f.findAt(((Length/2., -Width/2., 0.),),((Length/2., -h/2., 0.),))
    f_bottom = f.findAt(((Length/2., -Width, Depth/2.),))
    f_side1 = f.findAt(((0, -Width/2., Depth/2.),), ((0, -h/2., Depth/2.),))
    f_side2 = f.findAt(((Length, -Width/2., Depth/2.),), ((Length, -h/2., Depth/2.),))
    f_cortex_top = f.findAt(((Length/2., 0.0, Depth/2),))
    # f_cortex_bottom = f.findAt(((53.333333, -h, 0.083333), ))
    f_interface = f.findAt(((Length/2., -h, Depth/2),))

    # Define coordinate of sets-edges
    e_top = e.findAt(((Length/2., 0.0, Depth), ))  
    e_interface = e.findAt(((Length/2., -h, Depth),), ((Length/2., -h, 0.0),), )   
    e_interface_2 = e.findAt(((Length/2., -h, Depth), ))

    # Define coordinate of the sets-surfaces
    s_top = f.findAt(((Length/2., 0.0, Depth/2),))
    s_cortex = f.findAt(((Length/2., -h/2., Depth),))
    s_subcortex = f.findAt(((Length/2., -(Width-h)/2., Depth),))
    s_side1 = f.findAt(((0, -Width/2., Depth/2.),), ((0, -h/2., Depth/2.),))
    s_side2 = f.findAt(((Length, -Width/2., Depth/2.),), ((Length, -h/2., Depth/2.),))
    s_interface = f.findAt(((Length/2., -h, Depth/2),))


    # Define coordinate of sets-cells
    cells = c.findAt(((Length/2., -Width/2., Depth/2),),((Length/2., -h/2., 0.),))

    # Assign sets-faces
    region = p.Set(faces=f_back, name='back') 
    region = p.Set(faces=f_bottom, name='bottom')
    region = p.Set(faces=f_side1, name='side-1')
    region = p.Set(faces=f_side2, name='side-2')
    region = p.Set(faces=s_cortex, name='cortex_front')
    region = p.Set(faces=s_subcortex, name='subcortex_front')
    region = p.Set(faces=f_cortex_top, name='cortex_top') 
    # region = p.Set(faces=f_cortex_bottom, name='cortex_bottom')
    region = p.Set(faces=f_interface, name='interface')
 
    # Assign sets-edges
    region = p.Set(edges=e_top, name='top-edge')
    region = p.Set(edges=e_interface, name='interface') 
    region = p.Set(edges=e_interface_2, name='interface-edge')

    # Assign sets-surfaces
    region = p.Surface(side1Faces=s_top, name='top_s')
    region = p.Surface(side1Faces=s_cortex, name='cortex_s')
    region = p.Surface(side1Faces=s_subcortex, name='subcortex_s')
    region = p.Surface(side1Faces=s_side1, name='side1_s')
    region = p.Surface(side1Faces=s_side2, name='side2_s')
    region = p.Surface(side1Faces=s_interface, name='interface_s')

    # Assign sets-cells
    region = p.Set(cells=cells, name='whole-domain')

    #Assign sets-node
    p.Set(vertices=verts, name='node')

    # Rigid body surface
    p1 = mdb.models[ModelName].parts['wall']

    r = p1.referencePoints
    refPoints=(r[2], )
    p1.Set(referencePoints=refPoints, name='ref')

#######################################################################################
#######################################################################################
def Create_Boundary_Conditions(ModelName, InstanceName, Step):
    a = mdb.models[ModelName].rootAssembly   

    # Add displacement and fixed bcs

    region = a.instances[InstanceName].sets['back']
    mdb.models[ModelName].DisplacementBC(name='back', createStepName='pressure', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    region = a.instances[InstanceName].sets['bottom']
    mdb.models[ModelName].DisplacementBC(name='bottom', createStepName='pressure', region=region, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    # mdb.models[ModelName].EncastreBC(name='bottom', createStepName=StepName, region=region, localCsys=None)

    
    region = a.instances[InstanceName].sets['side-1']
    mdb.models[ModelName].DisplacementBC(name='side-1', createStepName='pressure', region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    region = a.instances[InstanceName].sets['side-2']
    mdb.models[ModelName].DisplacementBC(name='side-2', createStepName='pressure', region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    region = a.instances[InstanceName].sets['cortex_front']
    mdb.models[ModelName].DisplacementBC(name='cortex_front', createStepName='pressure', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    region = a.instances[InstanceName].sets['subcortex_front']
    mdb.models[ModelName].DisplacementBC(name='subcortex_front', createStepName='pressure', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.instances['wall-1'].sets['ref']
    mdb.models[ModelName].EncastreBC(name='ref-left', createStepName=StepName, region=region, localCsys=None)

    region = a.instances['wall-2'].sets['ref']
    mdb.models[ModelName].EncastreBC(name='ref-right', createStepName=StepName, region=region, localCsys=None)

#######################################################################################
#######################################################################################     
def Create_Pressure(ModelName, InstanceName, Step, top_p):
    # Add surface traction
    a = mdb.models[ModelName].rootAssembly 

    v1 = a.instances[InstanceName].vertices
    e1 = a.instances[InstanceName].edges
    region = a.instances[InstanceName].surfaces['top_s']
    mdb.models[ModelName].SurfaceTraction(name='top-p', createStepName=StepName, region=region, magnitude=top_p, amplitude='Amp-1', directionVector=(v1[10], a.instances[InstanceName].InterestingPoint(edge=e1[15], rule=MIDDLE)), distributionType=UNIFORM, field='', localCsys=None, traction=GENERAL, resultant=ON)

#######################################################################################
####################################################################################### 
def Create_VP(ModelName, InstanceName, Step, VP):

    a = mdb.models[ModelName].rootAssembly 
    region = a.instances[InstanceName].surfaces['top_s']
    mdb.models[ModelName].Pressure(name='VP', createStepName=StepName, region=region, distributionType=VISCOUS, field='', magnitude=VP, amplitude='Amp-1')

#######################################################################################
#######################################################################################    
def Create_Contact(ModelName, InstanceName, Step, Dimensions):
    a = mdb.models[ModelName].rootAssembly

    mdb.models[ModelName].ContactProperty('IntProp-1') 
    mdb.models[ModelName].interactionProperties['IntProp-1'].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((10.0, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)
    mdb.models[ModelName].interactionProperties['IntProp-1'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
    mdb.models[ModelName].interactionProperties['IntProp-1'].GeometricProperties(contactArea=1.0, padThickness=0.1)

    mdb.models[ModelName].ContactProperty('IntProp-2')
    mdb.models[ModelName].interactionProperties['IntProp-2'].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((10.0, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)
    mdb.models[ModelName].interactionProperties['IntProp-2'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
    
    region = a.instances[InstanceName].surfaces['top_s']
    mdb.models[ModelName].SelfContactExp(name='Int-1', createStepName=StepName, surface=region, mechanicalConstraint=KINEMATIC, interactionProperty='IntProp-1')

    s1 = a.instances['wall-2'].faces
    side1Faces1 = s1.findAt(((0.0, -Width/2, 0.458333), ))
    region1=a.Surface(side1Faces=side1Faces1, name='left-wall')
    region2=a.instances[ModelName].surfaces['top_s']
    mdb.models[ ModelName].SurfaceToSurfaceContactExp(name ='Int-2', createStepName='Initial', main = region1, secondary = region2, mechanicalConstraint=KINEMATIC, sliding=FINITE, interactionProperty='IntProp-2', initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    s1 = a.instances['wall-1'].faces
    side1Faces1 = s1.findAt(((Length, -Width/2, 0.458333), ))
    region1=a.Surface(side2Faces=side1Faces1, name='right_wall')
    region2=a.instances[ModelName].surfaces['top_s']
    mdb.models[ModelName].SurfaceToSurfaceContactExp(name ='Int-3', createStepName='Initial', main = region1, secondary = region2, mechanicalConstraint=KINEMATIC, sliding=FINITE, interactionProperty='IntProp-2', initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
#######################################################################################    
#######################################################################################    
def Create_Mesh(ModelName, PartName, InstanceName, Step, Dimensions):  
    # Rigif plane mesh
    p = mdb.models[ModelName].parts['wall']
    p.seedPart(size=1.5, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models[ModelName].parts['wall']
    p.generateMesh()

    # Bilayer Mesh
    p = mdb.models[ModelName].parts[PartName]
    a = mdb.models[ModelName].rootAssembly   
    e = p.edges

    # Define meshing edges 
    e_subcortex1 = e.findAt(((Length, -Width/2., Depth  ), ), ((0.0, -Width/2., 0.0), ))
    e_subcortex2 = e.findAt(((Length, -Width/2., 0.0), ), ((0.0, -Width/2., Depth  ), ))
    e_cortex = e.findAt(((0.0, -h/2, 0.0), ), ((Length, -h/2, 0.0), ), (( 0.0, -h/2, Depth), ), ((Length, -h/2, Depth), ))
    e_width = e.findAt(((Length/2., -h, Depth), ), ((Length/2., -h, 0.0), ), ((Length/2., 0.0, 0.0), ), ((Length/2., -Width, Depth), ), ((Length/2, -Width, 0.0), ), ((Length/2, 0.0, Depth), ))
    
    # Assign number of elements to the edges
    p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=e_subcortex1, end2Edges=e_subcortex2, ratio=bias, number=esubcortex, constraint=FINER)
    p.seedEdgeByNumber(edges=e_cortex, number=ecortex, constraint=FINER)
    p.seedEdgeByNumber(edges=e_width, number=ewidth, constraint=FINER)

    # Define meshing edges 
    e_subcortex1 = e.findAt(((Length, -Width/2., Depth  ), ), ((0.0, -Width/2., 0.0), ))
    e_subcortex2 = e.findAt(((Length, -Width/2., 0.0), ), ((0.0, -Width/2., Depth  ), ))
    e_cortex = e.findAt(((0.0, -h/2, 0.0), ), ((Length, -h/2, 0.0), ), (( 0.0, -h/2, Depth), ), ((Length, -h/2, Depth), ))
    e_width = e.findAt(((Length/2., -h, Depth), ), ((Length/2., -h, 0.0), ), ((Length/2., 0.0, 0.0), ), ((Length/2., -Width, Depth), ), ((Length/2, -Width, 0.0), ), ((Length/2, 0.0, Depth), ))
    
    # Assign number of elements to the edges
    p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=e_subcortex1, end2Edges=e_subcortex2, ratio=bias, number=esubcortex, constraint=FINER)
    p.seedEdgeByNumber(edges=e_cortex, number=ecortex, constraint=FINER)
    p.seedEdgeByNumber(edges=e_width, number=ewidth, constraint=FINER)

    # Define element type & Create mesh
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl=ENHANCED, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    #elemType1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    #elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    #elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    
    c = p.cells
    c_all = c.findAt(((Length/2., -h/2., Depth), ), ((Length/2., -Width/2., Depth), ))
    pickedRegions =(c_all, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
    partInstances =(a.instances[InstanceName], )
    p.generateMesh(regions=partInstances)

    # Define perturbation
    nodes = p.sets['interface'].nodes
    coords = []
           
    for n in nodes:
        x = n.coordinates[0]
        y = n.coordinates[1]
        z = n.coordinates[2]

        percent = 0.02
        dlt = abs(x - Length/2)
        wavelength = Length/waves
        
	if ptype == 0:
		pv = 0.
	elif ptype == 1:
	    if dlt <= percent*Length/2.:
		pv = random.random() - 1.
	    else: 
		pv = 0.
	elif ptype == 2:
	    if dlt <= percent*Length/2.:
		pv = -cos(x/wavelength*2.*math.pi) 
	    else: 
		pv = 0.
	elif ptype == 3:
	    pv = 0
            wavelengths = []
	    if dlt <= percent*Length/2.:
                for i in range(waves):
                    wavelengths.append(random.random()*waves)   
		    pv = pv + sin(x/wavelengths[i]*2.*math.pi)
	    else: 
		pv = 0.
	y_p = y + perturb * pv
	coords.append((x, y_p, z))
    p.editNode(nodes=nodes, coordinates=coords) 

#######################################################################################
#######################################################################################
def Create_Output(ModelName, PartName, InstanceName, Step): 
    p = mdb.models[ModelName].parts[PartName]

    # Assign desired output to sets 
    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['whole-domain']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-2', createStepName=StepName, variables=('S', 'PE', 'LE', 'U', 'ELSE', 'COORD','EVOL','SDV'),  numIntervals=500, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['Cortex']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-3', createStepName=StepName, variables=('ELSE', 'COORD', 'U','EVOL','SDV'), numIntervals=500, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['Subcortex']
    mdb.models[ModelName].FieldOutputRequest(name='F-Output-4', createStepName=StepName, variables=('ELSE', 'COORD', 'U','EVOL','SDV'), numIntervals=500, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
        
    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['whole-domain']
    mdb.models[ModelName].HistoryOutputRequest(name='H-Output-1', createStepName='pressure', variables=('ALLIE', 'ALLKE'), numIntervals=1000, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

    regionDef=mdb.models[ModelName].rootAssembly.allInstances[InstanceName].sets['Cortex']
    mdb.models[ModelName].HistoryOutputRequest(name='H-Output-2', createStepName='pressure', variables=('ALLIE', 'ALLKE'), numIntervals=1000, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
#######################################################################################
#######################################################################################
def Create_Job(ModelName, JobName, USERSUB):

    mdb.Job(name=JobName, model=ModelName, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=USERSUB, scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=8, activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=8)
    mdb.jobs[JobName].writeInput(consistencyChecking=OFF)
    try: 
        mdb.jobs[JobName].submit(consistencyChecking=OFF)
        mdb.jobs[JobName].waitForCompletion()
    except: 
        print "Job %s crashed" % JobName

#######################################################################################
#######################################################################################
## for zero-pressure case j range should be 1 
for j in range(5):
    for i in range(31):

            # ======================================================
            # Dimensions
            # ======================================================
            Width = 20
            Length = 100
            Depth = 0.25
            partition = 0.05                                                # cortical thickness relative to whole thickness
            h = Width*partition                                             # cortical thickness
            perturb = 0.02*h
            waves = 1                                                       # one sinusoidal wave
            ptype = 2                                                       # 0 = no perturbation, 1 = random perturbation 2 = middle wave, 3 = superposition of waves
            Dimensions = [partition,h,Depth,Length,Width,waves,ptype,perturb] 
            # ======================================================
            # Mesh Parameters
            # ======================================================
            bias = 3                                                        # bias through width of subcortex
            ecortex = 4                                                     # number of elements in cortex 
            esubcortex = 24                                                 # number of elements in subcortex 
            edepth = 1                                                      # number of elements in depth
            ewidth = 240                                                    # number of elements in length
            Mesh = [bias, ecortex, esubcortex, edepth, ewidth]
            # ======================================================
            # Step Parameters
            # ======================================================
            StepName = 'pressure'
            TotalTime = 1
            MS = 100 # Dont change this one
            Step = [StepName, TotalTime, MS]

            # ======================================================
            # Material Parameters
            # ======================================================
            CorticalMaterial = 'GRAY'
            SubCortMaterial = 'WHIT'
            mu_subcortex = 1e-3
            bulk_subcortex = 100*mu_subcortex # Dont change this one      
            lambda_subcortex = bulk_subcortex-(2*mu_subcortex/3) 
            Gctx =  2.0 
            density = 1e-11 # Dont change this one
            ratio = [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60, 70, 80 ,90, 100, 200, 400, 600, 800, 1000]
            Materials = [CorticalMaterial, SubCortMaterial, mu_subcortex, lambda_subcortex, Gctx, density, ratio]

            # ======================================================
            # Boundary Conditions 
            # ======================================================
            pressure_ratio = [0.5*1e-3, 1.0*1e-3, 2.0*1e-3, 3.0*1e-3, 4.0*1e-3]
            top_p = ratio[i]*pressure_ratio[j]
            VP = 1e-5 # Dont change this one
            # ======================================================
                # Model
            # ======================================================
            ModelName = 'bilayer-%d' %(i)
            PartName = 'bilayer-%d' %(i)
            InstanceName = 'bilayer-%d' %(i)
            JobName = 'p%d-%d' %(j,i)
            # ======================================================
                # User Material
            # ======================================================
            UMAT = "../JAM23_growth_VUMAT.f" # Change this to your own directory 
            USERSUB = UMAT #= os.getcwd()
            # ======================================================
                # Call Functions
            # ======================================================
            Create_Bilayered_Rectangle(ModelName, PartName, Dimensions)
            Create_Rigid_wall(ModelName)
            Create_Step(ModelName, Step)
            Create_Material(ModelName, Materials, ratio[i])
            Create_Section(ModelName, PartName, Dimensions, Materials)
            Create_Assembly(ModelName, PartName, InstanceName, Dimensions)
            Create_Sets(ModelName, PartName, Dimensions)
            Create_Boundary_Conditions(ModelName, InstanceName, Step)
            # Create_Pressure(ModelName, InstanceName, Step, top_p) ## Comment out for non zero-pressure cases 
            Create_VP(ModelName, InstanceName, Step, VP)
            Create_Contact(ModelName, InstanceName, Step, Dimensions)   
            Create_Mesh(ModelName, PartName, InstanceName, Step, Dimensions)
            Create_Output(ModelName, PartName, InstanceName, Step)
            Create_Job(ModelName, JobName, USERSUB)
