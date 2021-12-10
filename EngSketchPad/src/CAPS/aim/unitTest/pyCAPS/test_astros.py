from __future__ import print_function
import unittest

import os
import glob
import shutil

import sys

import pyCAPS

try:
    ASTROS_ROOT = os.environ["ASTROS_ROOT"]
except KeyError:
   print("Please set the environment variable ASTROS_ROOT")

# Helper function to check if an executable exists
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    try:
        exe_file = os.path.join(ASTROS_ROOT, "astros.exe")
        if is_exe(exe_file):
            return exe_file
    except:
        pass

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


class TestAstros(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.problemName = "workDir_astros"
        cls.iProb = 1
        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)
            
    # This is executed prior to each test
    def setUp(self):
        self.astros_exe = which("astros.exe")
        if self.astros_exe == None:
            self.skipTest("No astros.exe executable")
            return

    def run_astros(self, astros):

        Proj_Name = astros.input.Proj_Name

        # Run AIM pre-analysis
        astros.preAnalysis()

        ####### Run astros####################
        print ("\n\nRunning Astros......")

        # Copy files needed to run astros
        astros_files = ["ASTRO.D01",  # *.DO1 file
                        "ASTRO.IDX"]  # *.IDX file
        for file in astros_files:
            if not os.path.isfile(file):
                shutil.copy2(ASTROS_ROOT + os.sep + file, os.path.join(astros.analysisDir, file))

        # Run Astros via system call
        astros.system(self.astros_exe + " < " + Proj_Name +  ".dat > " + Proj_Name + ".out");

        print ("Done running Astros!")
        ########################################

        # Run AIM post-analysis
        astros.postAnalysis()


    def test_Plate(self):

        filename = os.path.join("..","csmData","feaSimplePlate.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=filename, outLevel=0); self.__class__.iProb += 1

        mesh = myProblem.analysis.create(aim="egadsTessAIM")
        
        astros = myProblem.analysis.create(aim = "astrosAIM",
                                           name = "astrosPlate")

        mesh.input.Edge_Point_Min = 3
        mesh.input.Edge_Point_Max = 4

        mesh.input.Mesh_Elements = "Quad"

        mesh.input.Tess_Params = [.25,.01,15]
        
        # Link the mesh
        astros.input["Mesh"].link(mesh.output["Surface_Mesh"])

        # Set analysis type
        astros.input.Analysis_Type = "Static"

        # Set materials
        madeupium    = {"materialType" : "isotropic",
                        "youngModulus" : 72.0E9 ,
                        "poissonRatio": 0.33,
                        "density" : 2.8E3}

        astros.input.Material = {"Madeupium": madeupium}

        # Set properties
        shell  = {"propertyType" : "Shell",
                  "membraneThickness" : 0.006,
                  "material"        : "madeupium",
                  "bendingInertiaRatio" : 1.0, # Default
                  "shearMembraneRatio"  : 5.0/6.0} # Default

        astros.input.Property = {"plate": shell}

        # Set constraints
        constraint = {"groupName" : "plateEdge",
                      "dofConstraint" : 123456}

        astros.input.Constraint = {"edgeConstraint": constraint}

        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : 2.e6}

        # Set loads
        astros.input.Load = {"appliedPressure": load }


        astros.input.File_Format = "Small"
        astros.input.Mesh_File_Format = "Small"
        astros.input.Proj_Name = "astrosPlateSmall"

        # Run Small format
        astros.runAnalysis()

        astros.input.File_Format = "Large"
        astros.input.Mesh_File_Format = "Large"
        astros.input.Proj_Name = "astrosPlateLarge"

        # Run Large format
        astros.runAnalysis()

        astros.input.File_Format = "Free"
        astros.input.Mesh_File_Format = "Free"
        astros.input.Proj_Name = "astrosPlateFree"

        # Run Free format
        astros.runAnalysis()

    def test_Aeroelastic(self):

        filename = os.path.join("..","csmData","feaWingOMLAero.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=filename, outLevel=0); self.__class__.iProb += 1

        mesh =  myProblem.analysis.create(aim="egadsTessAIM")
        
        astros =  myProblem.analysis.create(aim = "astrosAIM",
                                            name = "astrosAero",
                                            autoExec = False)

        astros.input.Proj_Name = "astrosAero"

        mesh.input.Edge_Point_Min = 3
        mesh.input.Edge_Point_Max = 4

        mesh.input.Mesh_Elements = "Quad"

        mesh.input.Tess_Params = [.25,.01,15]

        # Link the mesh
        astros.input["Mesh"].link(mesh.output["Surface_Mesh"])

        # Set analysis type
        astros.input.Analysis_Type = "Aeroelastic"

        # Set analysis
        trim1 = { "analysisType" : "AeroelasticStatic",
                  "trimSymmetry" : "SYM",
                  "analysisConstraint" : ["ribConstraint"],
                  "analysisSupport" : ["ribSupport"],
                  "machNumber"     : 0.5,
                  "dynamicPressure": 50000,
                  "density" : 1.0,
                  "rigidVariable"  : "ANGLEA",
                  "rigidConstraint": ["URDD3"],
                  "magRigidConstraint" : [-50],
                  }


        astros.input.Analysis = {"Trim1" : trim1}

        # Set materials
        unobtainium  = {"youngModulus" : 2.2E6 ,
                        "poissonRatio" : .5,
                        "density"      : 7850}

        astros.input.Material = {"Unobtainium" : unobtainium}

        # Set property
        shell  = {"propertyType"        : "Shell",
                  "membraneThickness"   : 0.2,
                  "bendingInertiaRatio" : 1.0, # Default
                  "shearMembraneRatio"  : 5.0/6.0} # Default }

        shell2  = {"propertyType"       : "Shell",
                  "membraneThickness"   : 0.002,
                  "bendingInertiaRatio" : 1.0, # Default
                  "shearMembraneRatio"  : 5.0/6.0} # Default }

        astros.input.Property = {"Rib_Root": shell,
                                 "Skin": shell2}


        # Defined Connections
        connection = {"dofDependent"   : 123456,
                      "connectionType" : "RigidBody"}

        astros.input.Connect = {"Rib_Root": connection}


        # Set constraints
        constraint = {"groupName"     : ["Rib_Root_Point"],
                      "dofConstraint" : 12456}

        astros.input.Constraint = {"ribConstraint": constraint}

        # Set supports
        support = {"groupName" : ["Rib_Root_Point"],
                   "dofSupport": 3}

        astros.input.Support = {"ribSupport": support}

        # Aero
        wing = {"groupName"    : "Wing",
                "numChord"     : 8,
                "numSpanTotal" : 24}

        # Note the surface name corresponds to the capsBound found in the *.csm file. This links
        # the spline for the aerodynamic surface to the structural model
        astros.input.VLM_Surface = {"Skin_Top": wing}

        # Run
        self.run_astros(astros)


if __name__ == '__main__':
    unittest.main()
