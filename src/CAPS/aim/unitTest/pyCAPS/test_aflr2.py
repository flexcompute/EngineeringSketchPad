from __future__ import print_function
import unittest

import os
import glob
import shutil

import pyCAPS

class Testaflr2(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_aflr2Test"
        cls.iDir = 1
        cls.iProb = 1
        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove problem directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

    def test_setInput(self):

        file = "../csmData/cfd2D.csm"
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr2AIM")

        myAnalysis.input.Proj_Name = "pyCAPS_aflr2_Test"
        myAnalysis.input.Tess_Params = [.25,.01,15]
        myAnalysis.input.Mesh_Quiet_Flag = True
        myAnalysis.input.Mesh_Format = "Tecplot"
        myAnalysis.input.Mesh_ASCII_Flag = True
        myAnalysis.input.Mesh_Gen_Input_String = "mquad=1 mpp=3"
        myAnalysis.input.Edge_Point_Min = 5
        myAnalysis.input.Edge_Point_Max = 300

        airfoil = {"numEdgePoints" : 100, "edgeDistribution" : "Tanh", "initialNodeSpacing" : [0.001, 0.001]}
        myAnalysis.input.Mesh_Sizing = {"Airfoil"   : airfoil,
                                        "TunnelWall": {"numEdgePoints" : 50},
                                        "InFlow"    : {"numEdgePoints" : 25},
                                        "OutFlow"   : {"numEdgePoints" : 25},
                                        "2DSlice"   : {"tessParams" : [0.50, .01, 45]}}

        # Run
        myAnalysis.runAnalysis()

    def test_reenter(self):

        file = "../csmData/cfd2D.csm"
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr2AIM")

        myAnalysis.input.Mesh_Quiet_Flag = True

        myAnalysis.input.Mesh_Format = "Tecplot"

        airfoil = {"numEdgePoints" : 100, "edgeDistribution" : "Tanh", "initialNodeSpacing" : [0.001, 0.001]}
        myAnalysis.input.Mesh_Sizing = {"Airfoil"   : airfoil,
                                        "TunnelWall": {"numEdgePoints" : 50},
                                        "InFlow"    : {"numEdgePoints" : 25},
                                        "OutFlow"   : {"numEdgePoints" : 25},
                                        "2DSlice"   : {"tessParams" : [0.50, .01, 45]}}

        myAnalysis.input.Proj_Name = "pyCAPS_aflr2_Tri"

        # Run 1st time to make triangle mesh
        myAnalysis.runAnalysis()

        myAnalysis.input.Mesh_Gen_Input_String = "mquad=1 mpp=3"
        myAnalysis.input.Proj_Name = "pyCAPS_aflr2_Quad"

        # Run 2nd time to make quad mesh
        myAnalysis.runAnalysis()

if __name__ == '__main__':
    unittest.main()
