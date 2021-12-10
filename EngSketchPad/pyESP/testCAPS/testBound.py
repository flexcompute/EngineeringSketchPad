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
        cls.analysisDir = "UnitTest"
        cls.projectName = "basicTest"
        cls.iProb = 1

        cls.cleanUp()

        cls.boundName = "test"

        cls.myProblem = pyCAPS.capsProblem()

        cls.myGeometry = cls.myProblem.loadCAPS(cls.file, cls.projectName, verbosity=0)

        cls.myAnalysis1 = cls.myProblem.loadAIM(aim = "fun3dAIM",
                                               analysisDir = cls.analysisDir + "1")
        
        cls.myAnalysis2 = cls.myProblem.loadAIM(aim = "nastranAIM",
                                               analysisDir = cls.analysisDir + "2")

    @classmethod
    def tearDownClass(cls):
        del cls.myGeometry
        del cls.myAnalysis1
        del cls.myAnalysis2
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.projectName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        # Remove default projectName
        base = os.path.basename(__main__.__file__)
        projectName = os.path.splitext(base)[0]
        if os.path.isdir(projectName):
            shutil.rmtree(projectName)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

#         if os.path.isfile("myAnalysisGeometry.egads"):
#             os.remove("myAnalysisGeometry.egads")

         # \param capsBound Name of capsBound to use for the data bound. 
    # \param variableName Single or list of variables names to add. 
    # \param aimSrc Single or list of AIM names that will be the data sources for the bound.
    # \param aimDest Single or list of AIM names that will be the data destinations during the transfer. 
    # \param transferMethod Single or list of transfer methods to use during the transfer.
    # \param initValueDest Single or list of initial values for the destication data.
    #
    # \return Optionally returns the reference to the data bound dictionary (\ref dataBound) entry created 
    # for the bound class object (\ref pyCAPS.capsBound).
    
    # Test bound creation
    def test_boundInit(self):
        name = self.boundName + "1"

        myBound= self.myProblem.createDataBound(capsBound = name,
                                                variableName = "Pressure",
                                                aimSrc = self.myAnalysis1.aimName,
                                                aimDest = self.myAnalysis2.aimName)

        self.assertEqual(myBound, self.myProblem.dataBound[name])
        
        # Repeat bound name
        with self.assertRaises(pyCAPS.CAPSError) as e:
             myBound= self.myProblem.createDataBound(capsBound = name,
                                                     variableName = "Pressure",
                                                     aimSrc = self.myAnalysis1.aimName,
                                                     aimDest = self.myAnalysis2.aimName)
             
        
    
    # Test bound creation
    def test_boundInit2(self):
        name = self.boundName + "2"
        
        myBound = pyCAPS.capsBound(self.myProblem, 
                                   capsBound = name,
                                   variableName = "Pressure",
                                   aimSrc = self.myAnalysis1.aimName,
                                   aimDest = self.myAnalysis2.aimName)

        self.assertEqual(myBound, self.myProblem.dataBound[name])
        
    # Test bound creation with analysis object
    def test_boundInit2(self):
        name = self.boundName + "3"
        
        myBound = pyCAPS.capsBound(self.myProblem, 
                                   capsBound = name,
                                   variableName = "Pressure",
                                   aimSrc = self.myAnalysis1,
                                   aimDest = self.myAnalysis2.aimName)

        self.assertEqual(myBound, self.myProblem.dataBound[name])
    
    # Test bound creation without a problem object
    def test_boundInit3(self):
     
        name = self.boundName + "Dummy"

        with self.assertRaises(TypeError) as e:
            myBound = pyCAPS.capsBound("Not an problme", 
                                       capsBound = name,
                                       variableName = "Pressure",
                                       aimSrc = self.myAnalysis1.aimName,
                                       aimDest = self.myAnalysis2.aimName)

        
    # Test lack of inputs 
    def test_noInputs(self):
        name = self.boundName + "Dummy"

        # no name
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myBound= self.myProblem.createDataBound(variableName = "Pressure",
                                                    aimSrc = self.myAnalysis1.aimName,
                                                    aimDest = self.myAnalysis2.aimName)

        # variable name 
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myBound= self.myProblem.createDataBound(capsBound = name,
                                                    aimSrc = self.myAnalysis1.aimName,
                                                    aimDest = self.myAnalysis2.aimName)

        # aim source name 
  
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myBound= self.myProblem.createDataBound(capsBound = name,
                                                    variableName = "Pressure",
                                                    aimDest = self.myAnalysis2.aimName)

    # Test mismatch of inputs 
    def test_mismatchInputs(self):
        name = self.boundName + "Dummy"

        # Src
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myBound = pyCAPS.capsBound(self.myProblem, 
                                       capsBound = name,
                                       variableName = ["Pressure",  "Displacement"],
                                       aimSrc = self.myAnalysis1.aimName,
                                       aimDest = self.myAnalysis2.aimName)

        # Dest
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myBound = pyCAPS.capsBound(self.myProblem, 
                                       capsBound = name,
                                       variableName = "Pressure",
                                       aimSrc = [self.myAnalysis1.aimName, self.myAnalysis2.aimName],
                                       aimDest = self.myAnalysis2.aimName)
    
    # Bad analysis name
    def test_invalidAnalysisName(self):
        name = self.boundName + "Dummy"

        # aim source name 
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myBound = pyCAPS.capsBound(self.myProblem, 
                                          capsBound = name,
                                          variableName = "Pressure",
                                          aimSrc = "x",
                                          aimDest = self.myAnalysis2.aimName)
            
        # aim destination name 
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myBound = pyCAPS.capsBound(self.myProblem, 
                                          capsBound = name,
                                          variableName = "Pressure",
                                          aimSrc = self.myAnalysis1.aimName,
                                          aimDest = "Y")
            

    # Test bound to make sure variables in bound are being set correctly 
    def test_variableSet(self):
        name = self.boundName + "4"
        
        myBound = pyCAPS.capsBound(self.myProblem, 
                                   capsBound = name,
                                   variableName = "Pressure",
                                   aimSrc = self.myAnalysis1.aimName,
                                   aimDest = self.myAnalysis2.aimName)

        self.assertEqual(myBound.variables, ["Pressure"])
        
        if pyVersion.major > 2:
            self.assertCountEqual(list(myBound.vertexSet.keys()), [self.myAnalysis1.aimName, self.myAnalysis2.aimName]) 
        else: 
            self.assertItemsEqual(list(myBound.vertexSet.keys()), [self.myAnalysis1.aimName, self.myAnalysis2.aimName]) 

if __name__ == '__main__':
    unittest.main() 
