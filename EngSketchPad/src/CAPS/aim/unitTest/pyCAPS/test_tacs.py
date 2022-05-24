import unittest

import os
import glob
import shutil

import sys

import pyCAPS

class TestTACS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.problemName = "workDir_tacs"
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

    def test_PlateSensFile(self):

        filename = os.path.join("..","csmData","feaSimplePlate.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=filename, outLevel=0); self.__class__.iProb += 1

        mesh = myProblem.analysis.create(aim="egadsTessAIM")

        tacs = myProblem.analysis.create(aim = "tacsAIM",
                                         name = "tacs")

        mesh.input.Edge_Point_Min = 3
        mesh.input.Edge_Point_Max = 4

        mesh.input.Mesh_Elements = "Quad"

        mesh.input.Tess_Params = [.25,.01,15]

        NumberOfNode = mesh.output.NumberOfNode

        # Link the mesh
        tacs.input["Mesh"].link(mesh.output["Surface_Mesh"])

        # Set analysis type
        tacs.input.Analysis_Type = "Static"

        # Set materials
        madeupium    = {"materialType" : "isotropic",
                        "youngModulus" : 72.0E9 ,
                        "poissonRatio": 0.33,
                        "density" : 2.8E3}

        tacs.input.Material = {"Madeupium": madeupium}

        # Set properties
        shell  = {"propertyType" : "Shell",
                  "membraneThickness" : 0.006,
                  "material"        : "madeupium",
                  "bendingInertiaRatio" : 1.0, # Default
                  "shearMembraneRatio"  : 5.0/6.0} # Default

        tacs.input.Property = {"plate": shell}

        # Set constraints
        constraint = {"groupName" : "plateEdge",
                      "dofConstraint" : 123456}

        tacs.input.Constraint = {"edgeConstraint": constraint}

        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : 2.e6}

        # Set loads
        tacs.input.Load = {"appliedPressure": load }


        tacs.input.File_Format = "Small"
        tacs.input.Mesh_File_Format = "Small"
        tacs.input.Proj_Name = "astrosPlateSmall"

        # make a unit place so the snensitvity is simply the Cartesian coordinate
        myProblem.geometry.despmtr.plateLength = 1.0
        myProblem.geometry.despmtr.plateWidth = 1.0

        # Setup an AnalysisIn design variable
        thickness = 0.1;
        desvar = {"groupName" : "plate",
                  "initialValue" : thickness,
                  "lowerBound" : thickness*0.5,
                  "upperBound" : thickness*1.5,
                  "maxDelta"   : thickness*0.1}

        tacs.input.Design_Variable = {"plateLength" : {},
                                      "thick1" : desvar}
        
        desvarR = {"variableType": "Property",
                   "fieldName" : "T",
                   "constantCoeff" : 0.0,
                   "groupName" : "thick1",
                   "linearCoeff" : 1.0}

        tacs.input.Design_Variable_Relation = {"thick1" : desvarR}

        for inode in range(NumberOfNode):
            tacs.input.File_Format = "Small"

            # Write input files
            tacs.preAnalysis()

            if inode == 0:
                xyz = [None]*NumberOfNode
                filename = os.path.join(tacs.analysisDir, tacs.input.Proj_Name+".bdf")
                with open(filename, "r") as f:
                    f.readline()
                    for i in range(NumberOfNode):
                        line = f.readline()
                        line = line.split()
                        xyz[i] = [float(line[-3]), float(line[-2]), float(line[-1])]

            # Create a dummy sensitivity file
            filename = os.path.join(tacs.analysisDir, tacs.input.Proj_Name+".sens")
            with open(filename, "w") as f:
                f.write("2\n")                       # Two functionals
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

            tacs.postAnalysis()

            Func1             = tacs.dynout["Func1"].value
            Func1_plateLength = tacs.dynout["Func1"].deriv("plateLength")

            self.assertEqual(Func1, 42*(inode+1))
            self.assertAlmostEqual(Func1_plateLength, xyz[inode][0], 5)

            Func2             = tacs.dynout["Func2"].value
            Func2_plateLength = tacs.dynout["Func2"].deriv("plateLength")

            self.assertEqual(Func2, 21*(inode+1))
            self.assertAlmostEqual(Func2_plateLength, xyz[inode][0], 5)


        tacs.input.Design_Variable = {"plateWidth" : {}}
        tacs.input.Design_Variable_Relation = None

        for inode in range(NumberOfNode):
            tacs.input.File_Format = "Small"

            # Write input files
            tacs.preAnalysis()

            if inode == 0:
                xyz = [None]*NumberOfNode
                filename = os.path.join(tacs.analysisDir, tacs.input.Proj_Name+".bdf")
                with open(filename, "r") as f:
                    f.readline()
                    for i in range(NumberOfNode):
                        line = f.readline()
                        line = line.split()
                        xyz[i] = [float(line[-3]), float(line[-2]), float(line[-1])]

            # Create a dummy sensitivity file
            filename = os.path.join(tacs.analysisDir, tacs.input.Proj_Name+".sens")
            with open(filename, "w") as f:
                f.write("2\n")                       # Two functionals
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

            tacs.postAnalysis()

            Func1             = tacs.dynout["Func1"].value
            Func1_plateWidth  = tacs.dynout["Func1"].deriv("plateWidth")

            self.assertEqual(Func1, 42*(inode+1))
            self.assertAlmostEqual(Func1_plateWidth, xyz[inode][1], 5)

            Func2             = tacs.dynout["Func2"].value
            Func2_plateWidth  = tacs.dynout["Func2"].deriv("plateWidth")

            self.assertEqual(Func2, 21*(inode+1))
            self.assertAlmostEqual(Func2_plateWidth, xyz[inode][1], 5)

        # Rotate so the 2nd coordinate in the plate is the z-coordinate
        myProblem.geometry.despmtr.angle = 90.0

        for inode in range(NumberOfNode):
            tacs.input.File_Format = "Small"

            # Write input files
            tacs.preAnalysis()

            if inode == 0:
                xyz = [None]*NumberOfNode
                filename = os.path.join(tacs.analysisDir, tacs.input.Proj_Name+".bdf")
                with open(filename, "r") as f:
                    f.readline()
                    for i in range(NumberOfNode):
                        line = f.readline()
                        line = line.split()
                        xyz[i] = [float(line[-3]), float(line[-2]), float(line[-1])]

            # Create a dummy sensitivity file
            filename = os.path.join(tacs.analysisDir, tacs.input.Proj_Name+".sens")
            with open(filename, "w") as f:
                f.write("2\n")                       # Two functionals
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

            tacs.postAnalysis()

            Func1             = tacs.dynout["Func1"].value
            Func1_plateWidth  = tacs.dynout["Func1"].deriv("plateWidth")

            self.assertEqual(Func1, 42*(inode+1))
            self.assertAlmostEqual(Func1_plateWidth, xyz[inode][2], 5)

            Func2             = tacs.dynout["Func2"].value
            Func2_plateWidth  = tacs.dynout["Func2"].deriv("plateWidth")

            self.assertEqual(Func2, 21*(inode+1))
            self.assertAlmostEqual(Func2_plateWidth, xyz[inode][2], 5)

if __name__ == '__main__':
    unittest.main()
