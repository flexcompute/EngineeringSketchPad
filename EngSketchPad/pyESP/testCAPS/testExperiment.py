import unittest

import os
import glob
import shutil
import __main__

from sys import version_info as pyVersion

import pyCAPS

class TestProblem(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.problemName = "testProblem"
        cls.iDir = 1
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

        # Remove default problemName
        base = os.path.basename(__main__.__file__)
        problemName = os.path.splitext(base)[0]
        if os.path.isdir(problemName):
            shutil.rmtree(problemName)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

        if os.path.isfile("myProblem.html"):
            os.remove("myProblem.html")

#==============================================================================
    # Create a bound
    def test_bound(self):

        problem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        fun3d = problem.analysis.create(aim = "fun3dAIM", name="fun3d")

        nastran = problem.analysis.create(aim = "nastranAIM", name="nastran")

        Upper_Left = problem.bound.create("Upper_Left")
        self.assertIs(Upper_Left, problem.bound["Upper_Left"])

        fun3dVert   = Upper_Left.vertexSet.create(fun3d, "fun3d")
        nastranVert = Upper_Left.vertexSet.create(nastran, "nastran")
        self.assertIs(fun3dVert  , Upper_Left.vertexSet["fun3d"])
        self.assertIs(nastranVert, Upper_Left.vertexSet["nastran"])

        fun3dVert.dataSet.create("Pressure", pyCAPS.fType.FieldOut)
        nastranVert.dataSet.create("Pressure", pyCAPS.fType.FieldIn)

        fun3dVert.dataSet.create("Displacement", pyCAPS.fType.FieldIn)
        nastranVert.dataSet.create("Displacement", pyCAPS.fType.FieldOut)

        fun3dVert.dataSet.create("EigenVector_1", pyCAPS.fType.FieldIn)
        nastranVert.dataSet.create("EigenVector_1", pyCAPS.fType.FieldOut)

        nastranVert.dataSet["Pressure"].link(fun3dVert.dataSet["Pressure"])
        fun3dVert.dataSet["Displacement"].link(nastranVert.dataSet["Displacement"])
        fun3dVert.dataSet["EigenVector_1"].link(nastranVert.dataSet["EigenVector_1"])

        Upper_Left.complete()

#         bounds = ["upperWing", "lowerWing", "leftTip", "riteTip"]
#         for capsBound in bounds:
#             problem.bound.create(capsBound)
#
#             astrosVert = problem.bound[capsBound].vertexSet.create(astros)
#             su2Vert  = problem.bound[capsBound].vertexSet.create(su2)
#
#             su2Vert.dataSet["Displacement"].link(astrosVert.dataSet["Displacement"], "Interpolate", init=(0,0,0))
#             astrosVert.dataSet["Pressure"].link(su2Vert.dataSet["Pressure"], "Conserve")
#
#
#         bounds = ["upperWing", "lowerWing", "leftTip", "riteTip"]
#         for capsBound in bounds:
#             problem.bound.create(capsBound)
#
#             astrosVert = problem.bound[capsBound].vertexSet.create(astros)
#             su2Vert  = problem.bound[capsBound].vertexSet.create(su2)
#
#             astrosVert.dataSet.create("Displacement", "Analysis")
#             su2Vert.dataSet.create("Displacement", "Interpolate", init=(0,0,0))
#
#             astrosVert.dataSet.create("Pressure", "Conserve")
#             su2Vert.dataSet.create("Pressure", "Analysis")
#
#             problem.bound[capsBound].complete()

#==============================================================================
    # Add an AIM
    def test_deprecations(self):

        # Initialize capsProblem object
        myProblem = pyCAPS.problem.capsProblem()

        # Load CSM file
        myGeometry = myProblem.loadCAPS(self.file, self.problemName+str(self.iProb), verbosity=0); self.__class__.iProb += 1

        # Load masstran aim
        masstran = myProblem.loadAIM(aim = "masstranAIM", altName = "masstran", capsIntent = "OML")

        # Min/Max number of points on an edge for quadding
        masstran.setAnalysisVal("Edge_Point_Min", 2)
        masstran.setAnalysisVal("Edge_Point_Max", 20)

        # Set materials
        madeupium    = {"materialType" : "isotropic",
                        "density"      : 10}
        unobtainium  = {"materialType" : "isotropic",
                        "density"      : 20}

        masstran.setAnalysisVal("Material", [("madeupium", madeupium),("unobtainium",unobtainium)])

        # Set properties
        shell1 = {"propertyType"      : "Shell",
                  "membraneThickness" : 2.0,
                  "material"          : "madeupium"}

        shell2 = {"propertyType"      : "Shell",
                  "membraneThickness" : 3.0,
                  "material"          : "unobtainium"}

        # set the propery for this analysis
        masstran.setAnalysisVal("Property", [("Wing1", shell1), ("Wing2", shell2)])

        # Only work with plate
        myGeometry.setGeometryVal("taper", 0.6)

        # compute Mass properties
        masstran.preAnalysis()
        masstran.postAnalysis()

        # properties of 1 plate
        Area     = masstran.getAnalysisOutVal("Area")
        Mass     = masstran.getAnalysisOutVal("Mass")
        Centroid = masstran.getAnalysisOutVal("Centroid")
        CG       = masstran.getAnalysisOutVal("CG")
        I        = masstran.getAnalysisOutVal("I_Vector")
        II       = masstran.getAnalysisOutVal("I_Tensor")


if __name__ == '__main__':
    unittest.main()
