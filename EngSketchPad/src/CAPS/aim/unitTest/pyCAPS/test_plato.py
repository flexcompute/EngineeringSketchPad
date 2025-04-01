import unittest

import os
import glob
import shutil

import sys

import pyCAPS

# get the value of _ESP_ROOT
try:
    _ESP_ROOT = os.environ["ESP_ROOT"]
except:
    raise RuntimeError("ESP_ROOT must be set -- Please fix the environment...")

class TestPLATO(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_plato"
        cls.iProb = 1
        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):
        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

    # This is executed prior to each test
    def setUp(self):
        if not os.path.isfile(os.path.join(_ESP_ROOT,"lib","platoAIM.so")):
            self.skipTest("Plato AIM not available")

#==============================================================================
    def test_PlateSensFile(self):

        filename = os.path.join("..","csmData","feaSimplePlate.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=filename, outLevel=0); self.__class__.iProb += 1

        # make a unit place so the snensitvity is simply the Cartesian coordinate
        myProblem.geometry.despmtr.plateLength = 1.0
        myProblem.geometry.despmtr.plateWidth = 1.0

        mesh = myProblem.analysis.create(aim="egadsTessAIM")

        plato = myProblem.analysis.create(aim = "platoAIM",
                                         name = "plato")

        mesh.input.Edge_Point_Min = 3
        mesh.input.Edge_Point_Max = 4

        mesh.input.Mesh_Elements = "Mixed"


        mesh.input.Mesh_Format = "Tecplot"

        mesh.input.Tess_Params = [.25,.01,15]

        NumberOfNode = mesh.output.NumberOfNode

        # Link the mesh
        plato.input["Mesh"].link(mesh.output["Surface_Mesh"])

        # Use mesh morphing
        plato.input.Mesh_Morph = True

        # Read sensitvity file
        plato.input.Design_SensFile = True

        # Setup an GeometryIn design variable
        plato.input.Design_Variable = {"plateLength" : {}}
        
        for inode in range(NumberOfNode):
            # Read sensitvity file
            plato.input.Design_SensFile = True
        
            # Write input files
            plato.preAnalysis()
        
            if inode == 0:
                xyz = [None]*NumberOfNode
                for i in range(NumberOfNode):
                    xyz[i] = [None,None,None]
                filename = os.path.join(mesh.analysisDir, mesh.input.Proj_Name+".tec")
                with open(filename, "r") as f:
                    f.readline() # Skip VARIABLES
                    f.readline() # Skip ZONE
                    # Read BLOCK format
                    for d in range(3):
                        for i in range(NumberOfNode):
                            xyz[i][d] = float(f.readline())
        
            # Create a dummy sensitivity file
            filename = os.path.join(plato.analysisDir, plato.input.Proj_Name+".sens")
            with open(filename, "w") as f:
                f.write("2 0\n")                     # Two functionals and no analysis design variable
                f.write("Func1\n")                   # 1st functional
                f.write("%f\n"%(42*(inode+1)))       # Value of Func1
                f.write("{}\n".format(NumberOfNode)) # Number Of Nodes for this functional
                for i in range(NumberOfNode):        # node ID, d(Func1)/d(xyz)
                    f.write("{} {} {} {}\n".format(i+1, 1 if i == inode else 0, 0, 0))
        
                f.write("Func2\n")                   # 2nd functiona;
                f.write("%f\n"%(21*(inode+1)))       # Value of Func2
                f.write("{}\n".format(NumberOfNode)) # Number Of Nodes for this functional
                for i in range(NumberOfNode):        # node ID, d(Func2)/d(xyz)
                    f.write("{} {} {} {}\n".format(i+1, 1 if i == inode else 0, 0, 0))
        
            plato.postAnalysis()
        
            Func1             = plato.dynout["Func1"].value
            Func1_plateLength = plato.dynout["Func1"].deriv("plateLength")
        
            self.assertEqual(Func1, 42*(inode+1))
            self.assertAlmostEqual(Func1_plateLength, xyz[inode][0], 5)
        
            Func2             = plato.dynout["Func2"].value
            Func2_plateLength = plato.dynout["Func2"].deriv("plateLength")
        
            self.assertEqual(Func2, 21*(inode+1))
            self.assertAlmostEqual(Func2_plateLength, xyz[inode][0], 5)

        # unlink the mesh to test Mesh_Morphing
        plato.input["Mesh"].unlink()

        myProblem.geometry.despmtr.plateWidth = 1.0

        plato.input.Design_Variable = {"plateWidth" : {}}

        for inode in range(NumberOfNode):
            # Read sensitvity file
            plato.input.Design_SensFile = True

            # Write input files
            plato.preAnalysis()

            if inode == 0:
                xyz = [None]*NumberOfNode
                for i in range(NumberOfNode):
                    xyz[i] = [None,None,None]
                filename = os.path.join(mesh.analysisDir, mesh.input.Proj_Name+".tec")
                with open(filename, "r") as f:
                    f.readline() # Skip VARIABLES
                    f.readline() # Skip ZONE
                    # Read BLOCK format
                    for d in range(3):
                        for i in range(NumberOfNode):
                            xyz[i][d] = float(f.readline())
            
            # Create a dummy sensitivity file
            filename = os.path.join(plato.analysisDir, plato.input.Proj_Name+".sens")
            with open(filename, "w") as f:
                f.write("2 0\n")                     # Two functionals and one analysis design variable
                f.write("Func1\n")                   # 1st functional
                f.write("%f\n"%(42*(inode+1)))       # Value of Func1
                f.write("{}\n".format(NumberOfNode)) # Number Of Nodes for this functional
                for i in range(NumberOfNode): # d(Func1)/d(xyz)
                    f.write("{} {} {} {}\n".format(i+1, 0, 1 if i == inode else 0, 0))

                f.write("Func2\n")                   # 2nd functiona;
                f.write("%f\n"%(21*(inode+1)))       # Value of Func2
                f.write("{}\n".format(NumberOfNode)) # Number Of Nodes for this functional
                for i in range(NumberOfNode): # d(Func2)/d(xyz)
                    f.write("{} {} {} {}\n".format(i+1, 0, 1 if i == inode else 0, 0))

            plato.postAnalysis()

            Func1             = plato.dynout["Func1"].value
            Func1_plateWidth  = plato.dynout["Func1"].deriv("plateWidth")

            self.assertEqual(Func1, 42*(inode+1))
            self.assertAlmostEqual(Func1_plateWidth, xyz[inode][1], 5)

            Func2             = plato.dynout["Func2"].value
            Func2_plateWidth  = plato.dynout["Func2"].deriv("plateWidth")

            self.assertEqual(Func2, 21*(inode+1))
            self.assertAlmostEqual(Func2_plateWidth, xyz[inode][1], 5)

        # Rotate so the 2nd coordinate in the plate is the z-coordinate
        myProblem.geometry.despmtr.angle = 90.0
        mesh.runAnalysis()
        
        plato.input.Design_Variable = {"plateWidth" : {}}

        for inode in range(NumberOfNode):
            # Read sensitvity file
            plato.input.Design_SensFile = True

            # Write input files
            plato.preAnalysis()

            if inode == 0:
                xyz = [None]*NumberOfNode
                for i in range(NumberOfNode):
                    xyz[i] = [None,None,None]
                filename = os.path.join(mesh.analysisDir, mesh.input.Proj_Name+".tec")
                with open(filename, "r") as f:
                    f.readline() # Skip VARIABLES
                    f.readline() # Skip ZONE
                    # Read BLOCK format
                    for d in range(3):
                        for i in range(NumberOfNode):
                            xyz[i][d] = float(f.readline())

            # Create a dummy sensitivity file
            filename = os.path.join(plato.analysisDir, plato.input.Proj_Name+".sens")
            with open(filename, "w") as f:
                f.write("2 0\n")                     # Two functionals and one analysis design variable
                f.write("Func1\n")                   # 1st functional
                f.write("%f\n"%(42*(inode+1)))       # Value of Func1
                f.write("{}\n".format(NumberOfNode)) # Number Of Nodes for this functional
                for i in range(NumberOfNode): # d(Func1)/d(xyz)
                    f.write("{} {} {} {}\n".format(i+1, 0, 0, 1 if i == inode else 0))

                f.write("Func2\n")                   # 2nd functiona;
                f.write("%f\n"%(21*(inode+1)))       # Value of Func2
                f.write("{}\n".format(NumberOfNode)) # Number Of Nodes for this functional
                for i in range(NumberOfNode): # d(Func2)/d(xyz)
                    f.write("{} {} {} {}\n".format(i+1, 0, 0, 1 if i == inode else 0))

            plato.postAnalysis()

            Func1             = plato.dynout["Func1"].value
            Func1_plateWidth  = plato.dynout["Func1"].deriv("plateWidth")

            self.assertEqual(Func1, 42*(inode+1))
            self.assertAlmostEqual(Func1_plateWidth, xyz[inode][2], 5)

            Func2             = plato.dynout["Func2"].value
            Func2_plateWidth  = plato.dynout["Func2"].deriv("plateWidth")

            self.assertEqual(Func2, 21*(inode+1))
            self.assertAlmostEqual(Func2_plateWidth, xyz[inode][2], 5)


if __name__ == '__main__':
    unittest.main()
