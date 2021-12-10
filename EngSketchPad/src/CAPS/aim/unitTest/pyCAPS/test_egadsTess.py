from __future__ import print_function
import unittest

import os
import glob
import shutil

import pyCAPS

class TestEGADS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_egadsTest"
        cls.iProb = 1
        cls.cleanUp()

        # Initialize a global Problem object
        cornerFile = os.path.join("..","csmData","cornerGeom.csm")
        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cornerFile, outLevel=0)

    @classmethod
    def tearDownClass(cls):
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

    def test_invalid_Mesh_Lenght_Scale(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        # Set project name so a mesh file is generated
        myAnalysis.input.Proj_Name = "pyCAPS_egadsTess_Test"
        
        # Set output grid format since a project name is being supplied - Tecplot  file
        myAnalysis.input.Mesh_Format = "Tecplot"

        # Set output grid format since a project name is being supplied - Tecplot  file
        myAnalysis.input.Mesh_Length_Factor = -1

        with self.assertRaises(pyCAPS.CAPSError) as e:
            myAnalysis.runAnalysis()

        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

    def test_setInput(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        myAnalysis.input.Proj_Name = "pyCAPS_EGADS_Test"
        myAnalysis.input.Mesh_Length_Factor = 1.05
        myAnalysis.input.Tess_Params = [0.1, 0.01, 15.0]
        myAnalysis.input.Mesh_Format = "Tecplot"
        myAnalysis.input.Mesh_ASCII_Flag = True
        myAnalysis.input.Edge_Point_Min = 2
        myAnalysis.input.Edge_Point_Max = 10
        
        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield": {"tessParams" : [0, 0.2, 30]}}
        myAnalysis.input.Mesh_Sizing = Mesh_Sizing

        myAnalysis.input.Mesh_Elements = "Quad"
        myAnalysis.input.Multiple_Mesh = True
        myAnalysis.input.TFI_Templates = False

    def test_SingleBody_AnalysisOutVal(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield": {"tessParams" : [0.3*80, 0.2*80, 30]}}
        myAnalysis.input.Mesh_Sizing = Mesh_Sizing

        # Execution is automatic

        # Assert AnalysisOutVals
        self.assertTrue(myAnalysis.output.Done)

        numNodes = myAnalysis.output.NumberOfNode
        self.assertGreater(numNodes, 0)
        numElements = myAnalysis.output.NumberOfElement
        self.assertGreater(numElements, 0)

    def test_MultiBody(self):

        file = os.path.join("..","csmData","cfdMultiBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")
        
        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield" : {"tessParams" : [0.3*80, 0.2*80, 30]}}
        myAnalysis.input.Mesh_Sizing = Mesh_Sizing

        # Run
        myAnalysis.runAnalysis()

    def test_reenter(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield": {"tessParams" : [0.3*80, 0.2*80, 30]}}
        myAnalysis.input.Mesh_Sizing = Mesh_Sizing

        myAnalysis.input.Mesh_Length_Factor = 1

        # Run 1st time
        myAnalysis.runAnalysis()

        myAnalysis.input.Mesh_Length_Factor = 2

        # Run 2nd time coarser
        myAnalysis.runAnalysis()


    def test_box(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Box",
                                                   capsIntent = ["box", "farfield"])

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

    def test_cylinder(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Cylinder",
                                                   capsIntent = ["cylinder", "farfield"])

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

    def test_cone(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Cone",
                                                   capsIntent = ["cone", "farfield"])

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

    def test_torus(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Torus",
                                                   capsIntent = ["torus", "farfield"])

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.viewGeometry()

    def test_sphere(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Sphere",
                                                   capsIntent = ["sphere", "farfield"])

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.viewGeometry()

    def test_boxhole(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "BoxHole",
                                                   capsIntent = ["boxhole", "farfield"])

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.viewGeometry()

    def test_bullet(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Bullet",
                                                   capsIntent = ["bullet", "farfield"])

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.viewGeometry()

    def test_all(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "All",
                                                   capsIntent = ["box", "cylinder", "cone", "torus", "sphere", "boxhole", "bullet", "farfield"]) 

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.viewGeometry()


if __name__ == '__main__':
    unittest.main()
