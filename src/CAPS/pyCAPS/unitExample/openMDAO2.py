from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Import OpenMDAO module
from openmdao.api import Problem, Group, IndepVarComp

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. 
print("Loading file into our capsProblem")
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm")

directory = "OpenMDAOUnitTest"
intent = "ALL"
# Load egadsTess aim 
myMeshSurf = myProblem.loadAIM(aim = "egadsTessAIM", 
                               analysisDir = directory,
                               capsIntent = intent)

# Load Tetgen aim 
myMeshVol = myProblem.loadAIM(aim = "tetgenAIM", 
                              analysisDir = directory,
                              capsIntent = intent,
                              parents = ["egadsTessAIM"])


# Load Fun3d aim 
fun3d = myProblem.loadAIM(aim = "fun3dAIM",  
                          analysisDir = "directory",
                          capsIntent = intent, 
                          parents = ["tetgenAIM"])

# Create OpenMDAOComponet - external code 
fun3d.createOpenMDAOComponent(["Mach"], 
                              ["CDtot"], 
                               executeCommand = ["nodet_mpi"]) # for avl execution 

# Set project name
fun3d.setAnalysisVal("Proj_Name", "fun3dTetgenTest")

fun3d.setAnalysisVal("Mesh_ASCII_Flag", False)
  
# Set AoA number 
fun3d.setAnalysisVal("Alpha", 1.0) 

# Set Mach number 
fun3d.setAnalysisVal("Mach", 0.5901)

# Set equation type
fun3d.setAnalysisVal("Eqn","compressible")

# Set Viscous term      
fun3d.setAnalysisVal("Viscous", "inviscid")

# Set number of iterations
fun3d.setAnalysisVal("Num_Iter",10)

# Set CFL number schedule
fun3d.setAnalysisVal("CFL_Schedule",[0.5, 3.0])

# Set read restart option
fun3d.setAnalysisVal("Restart_Read","off")

# Set CFL number iteration schedule
fun3d.setAnalysisVal("CFL_Schedule_Iter", [1, 40])
 
# Set overwrite fun3d.nml if not linking to Python library
fun3d.setAnalysisVal("Overwrite_NML", True)

##[setFUN3D]

# Set boundary conditions
inviscidBC1 = {"bcType" : "Inviscid", "wallTemperature" : 1}
inviscidBC2 = {"bcType" : "Inviscid", "wallTemperature" : 1.2}
fun3d.setAnalysisVal("Boundary_Condition", [("Wing1", inviscidBC1),
                                               ("Wing2", inviscidBC2),
                                               ("Farfield","farfield")])

# Close our CAPS pro
openMDAOProblem = Problem()
openMDAOProblem.root = Group()

# Create and connect inputs
openMDAOProblem.root.add('param1', IndepVarComp('Mach', 3.0))
openMDAOProblem.root.add('myAnalysis', fun3d.openMDAOComponent)

openMDAOProblem.root.connect('param1.Mach', 'myAnalysis.Mach')

# Run the ExternalCode Component
openMDAOProblem.setup()
openMDAOProblem.run_once()

print("Closing our CAPS problem")
myProblem.closeCAPS()