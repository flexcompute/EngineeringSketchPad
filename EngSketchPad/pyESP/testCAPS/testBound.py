import unittest

import os
import glob
import shutil
import __main__

from sys import version_info as pyVersion
from sys import version_info

import pyCAPS

class TestBound(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.problemName = "basicTest"
        cls.iProb = 1

        cls.cleanUp()

        cls.capsBound = "Upper_Left"

        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file, outLevel=0)

        cls.su2 = cls.myProblem.analysis.create(aim = "su2AIM")
        
        cls.astros = cls.myProblem.analysis.create(aim = "astrosAIM")

    @classmethod
    def tearDownClass(cls):
        cls.myProblem.close()
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

#==============================================================================
    # Test bound creation
    def test_boundInit(self):
        
        myProblem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile=self.file, outLevel=0); self.iProb += 1

        su2    = myProblem.analysis.create(aim = "su2AIM")
        astros = myProblem.analysis.create(aim = "astrosAIM")

        bound = myProblem.bound.create(capsBound = self.capsBound)

        self.assertEqual(bound, myProblem.bound[self.capsBound])
        
        # Repeat bound name is not allowed
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myProblem.bound.create(capsBound = self.capsBound)
            
        # Create the vertex sets on the bound for su2 and astros analysis
        su2Vset    = bound.vertexSet.create(self.su2)
        astrosVset = bound.vertexSet.create(self.astros)

        self.assertEqual(su2Vset   , bound.vertexSet[self.su2.name])
        self.assertEqual(astrosVset, bound.vertexSet[self.astros.name])

        # Create pressure data sets
        su2_Pressure    = su2Vset.dataSet.create("Pressure", pyCAPS.fType.FieldOut)
        astros_Pressure = astrosVset.dataSet.create("Pressure", pyCAPS.fType.FieldIn)

        self.assertEqual(su2_Pressure   , su2Vset.dataSet["Pressure"])
        self.assertEqual(astros_Pressure, astrosVset.dataSet["Pressure"])

        # Create displacement data sets
        su2_Displacement    = su2Vset.dataSet.create("Displacement", pyCAPS.fType.FieldIn, init=[0,0,0])
        astros_Displacement = astrosVset.dataSet.create("Displacement", pyCAPS.fType.FieldOut)

        self.assertEqual(su2_Displacement   , su2Vset.dataSet["Displacement"])
        self.assertEqual(astros_Displacement, astrosVset.dataSet["Displacement"])

        # Link the data sets
        astros_Pressure.link(su2_Pressure, "Conserve")
        su2_Displacement.link(astros_Displacement, "Interpolate")

        # Close the bound as complete (cannot create more vertex or data sets)
        bound.close()
        
        # Explicitly close the problem to test the cleanup process
        myProblem.close()


if __name__ == '__main__':
    unittest.main() 
