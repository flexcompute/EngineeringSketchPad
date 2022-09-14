import unittest

import os
import glob
import math
import shutil
import __main__

from sys import version_info as pyVersion
from sys import version_info

from pyCAPS import caps

class TestAnalysis(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.fileWrite = "unitGeomAnalysis.egads"
        cls.analysisDir = "UnitTest"
        cls.projectName = "testProblem"
        cls.iDir = 1
        cls.iProb = 1
        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
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

        # Remove created files
        if os.path.isfile(cls.fileWrite):
            os.remove(cls.fileWrite)

#=============================================================================-
    # Same AIM twice
    def test_makeAnalysis(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        analysis1 = problem.makeAnalysis("fun3dAIM", name = self.analysisDir);
        analysis2 = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        with self.assertRaises(caps.CAPSError) as e:
            analysis3 = problem.makeAnalysis("fun3dAIM", name = self.analysisDir)

        self.assertEqual(e.exception.errorName, "CAPS_BADNAME")
 
#=============================================================================-
    # Duplicate AIM
#     def test_copyAIM(self):
#         problem = caps.open("TestProblem", self.file)
# 
#         analysis1 = problem.load("fun3dAIM", unitSys=None,
#                                   analysisDir = self.analysisDir)
#         analysis2 = analysis1.dupAnalysis()
# 
#         self.assertNotEqual(analysis2, None)

#=============================================================================-
    # Multiple intents
    def test_multipleIntents(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        analysis = None
        analysis = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir),
                                        intent = ["OML"]); self.__class__.iDir += 1

        self.assertNotEqual(analysis, None)
        self.assertEqual(2, analysis.size(caps.oType.BODIES, caps.sType.NONE))

        analysis = None
        analysis = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir),
                                        intent = ["CFD"]); self.__class__.iDir += 1

        self.assertNotEqual(analysis, None)
        self.assertEqual(3, analysis.size(caps.oType.BODIES, caps.sType.NONE))

#=============================================================================-
    # get all bodies in the analysis
    def test_getBodies(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        analysis = problem.makeAnalysis("fun3dAIM", name = self.analysisDir);

        bodies = analysis.getBodies()

        self.assertEqual(4, len(bodies))
 
#=============================================================================-
    # Set analysis value with units
    def test_setInputUnits(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1
        
        
        deg = caps.Unit("degree")
        rad = caps.Unit("radian")

        m    = caps.Unit("meter")
        kg   = caps.Unit("kg")
        s    = caps.Unit("s")
        K    = caps.Unit("Kelvin")

        unitSystem={"mass":kg, "length":m, "time":s, "temperature":K}
 
        analysis = problem.makeAnalysis("su2AIM", name = self.analysisDir + str(self.iDir), unitSys=unitSystem); self.__class__.iDir += 1
 
        betaObj = analysis.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, "Beta")
        betaObj.setValue(10.*math.pi/180. * rad)
        beta = betaObj.getValue()

        self.assertAlmostEqual(10, beta._value, 5) # Should be converted to degree
        self.assertEqual(caps.Unit("degree"), beta._units)

#=============================================================================-
    # Set/Unset analysis value
    def test_setUnsetInput(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1
 
        analysis = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        machObj = analysis.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, "Mach")
        
        data = machObj.getValue()
        self.assertEqual(None, data)
        
        machObj.setValue(0.5)
        data = machObj.getValue()
        self.assertEqual(0.5, data)

        machObj.setValue(None)
        data = machObj.getValue()
        self.assertEqual(None, data)

#=============================================================================-
    # Get a string and add it to another string -> Make sure byte<->str<->unicode is correct between Python 2 and 3
    def test_AnalysisString(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1
 
        analysis = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        Proj_NameObj = analysis.childByName(caps.oType.VALUE, caps.sType.ANALYSISIN, "Proj_Name")

        setString = "myTest"
        Proj_NameObj.setValue(setString)
 
        getString = Proj_NameObj.getValue()
        self.assertEqual(setString+"_TEST", getString+"_TEST")

#=============================================================================-
    # Get analysis out values
    def test_getAnalysisOutVal(self):
  
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1
 
        analysis = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        CLtotObj = analysis.childByName(caps.oType.VALUE, caps.sType.ANALYSISOUT, "CLtot")

        # Get all output values - analysis should be dirty since pre hasn't been run
        with self.assertRaises(caps.CAPSError) as e:
            data, units = CLtotObj.getValue()

        self.assertEqual(e.exception.errorName, "CAPS_DIRTY")

#=============================================================================-
    def test_parent1(self):
        
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        analysis1 = problem.makeAnalysis("egadsTessAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        analysis2 = problem.makeAnalysis("tetgenAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1
 
        analysis3 = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1
 
        self.assertNotEqual(None, analysis3)

#=============================================================================-
#     def test_parent2(self):
#         
#         name = "analysisInitTest_4"
#         index = "14"
#         
#         with self.assertRaises(TypeError) as e:
#                 
#             myAnalysis = pyCAPS.capsAnalysis(self.myProblem, 
#                                             aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + index)
#             
#             myAnalysis1 = pyCAPS.capsAnalysis(self.myProblem, 
#                                               aim = "egadsTessAIM",
#                                               analysisDir = self.analysisDir + index)
#                     
#             myAnalysis2 = pyCAPS.capsAnalysis(self.myProblem, 
#                                               aim = "tetgenAIM",
#                                               analysisDir = self.analysisDir + index, 
#                                               parents = [[], myAnalysis1], 
#                                               altName=name)
#                     
#             self.assertEqual(myAnalysis, self.myProblem.analysis[name])
#         
# #=============================================================================-
#     def test_parent3(self):
#         
#         name = "analysisInitTest_5"
#         index = "15"
#         
#         myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + index)
#         
#         myAnalysis1 = self.myProblem.loadAIM(aim = "egadsTessAIM",
#                                              analysisDir = self.analysisDir + index)
#                 
#         myAnalysis2 = self.myProblem.loadAIM(aim = "tetgenAIM",
#                                              analysisDir = self.analysisDir + index, 
#                                              parents = [myAnalysis.aimName, myAnalysis1], 
#                                              altName=name)
#                 
#         self.assertEqual(myAnalysis2, self.myProblem.analysis[name])
#      
# #=============================================================================-
#     def test_parent4(self):
#         
#         name = "analysisInitTest_6"
#         index = "16"
#         
#         with self.assertRaises(pyCAPS.CAPSError) as e:
#                 
#             myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                                 analysisDir = self.analysisDir + index)
#             
#             myAnalysis1 = self.myProblem.loadAIM(aim = "egadsTessAIM",
#                                                  analysisDir = self.analysisDir + index)
#                     
#             myAnalysis2 = self.myProblem.loadAIM(aim = "tetgenAIM",
#                                                  analysisDir = self.analysisDir + index, 
#                                                  parents = ["soda", myAnalysis1], 
#                                                  altName=name)
#                     
#             self.assertEqual(myAnalysis, self.myProblem.analysis[name])     
#     
# #=============================================================================-
#     def test_parent5(self):
#         
#         name = "analysisInitTest_7"
#         index = "17"
#         
#         myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + index)
#                 
#         myAnalysis2 = self.myProblem.loadAIM(aim = "tetgenAIM",
#                                              analysisDir = self.analysisDir + index, 
#                                              parents = myAnalysis, 
#                                              altName=name)
#                 
#         self.assertEqual(myAnalysis2, self.myProblem.analysis[name])  
 
#=============================================================================-
    # Get analysis info
    def test_getAnalysisInfo(self):
 
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        analysis = problem.makeAnalysis("fun3dAIM", 
                                        intent = ["CFD", "OML"],
                                        name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        dir, unitSys, major, minor, intent, fnames, franks, fInOut, execute, status = analysis.analysisInfo()
        self.assertTrue(os.path.isdir(dir))
        self.assertEqual(None, unitSys)

        self.assertEqual(["CFD", "OML"], intent)
        self.assertEqual(['Pressure', 'P', 'Cp', 'CoefficientOfPressure', 'Displacement', 'EigenVector', 'EigenVector_#'], fnames)
        self.assertEqual([1,1,1,1,3,3,3], franks)
        self.assertEqual([caps.fType.FieldOut,caps.fType.FieldOut,caps.fType.FieldOut,caps.fType.FieldOut,caps.fType.FieldIn,caps.fType.FieldIn,caps.fType.FieldIn], fInOut)
        self.assertEqual(0, execute)
        self.assertEqual(3, status)

#=============================================================================-
#     # Get analysis info with a NULL intent
#     def test_getAnalysisInfoIntent(self):
# 
#         myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + "5"
#                                             ) # Use default capsIntent
# 
#         infoDict = myAnalysis.getAnalysisInfo(printinfo = False,  infoDict = True)
#         self.assertEqual(infoDict["intent"], None)
# 
#         myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + "7",
#                                             capsIntent = None) # Use no capsIntent
# 
#         infoDict = myAnalysis.getAnalysisInfo(printinfo = False,  infoDict = True)
#         self.assertEqual(infoDict["intent"], None)
# 
#         myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + "8",
#                                             capsIntent = "") # Use no capsIntent
# 
#         infoDict = myAnalysis.getAnalysisInfo(printinfo = False,  infoDict = True)
#         self.assertEqual(infoDict["intent"], None)
# 
#         myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + "9",
#                                             capsIntent = "STRUCTURE") # Use different capsIntent
# 
#         infoDict = myAnalysis.getAnalysisInfo(printinfo = False,  infoDict = True)
#         self.assertEqual(infoDict["intent"], "STRUCTURE")
# 
#         myAnalysis = self.myProblem.loadAIM(aim = "fun3dAIM",
#                                             analysisDir = self.analysisDir + "10",
#                                             capsIntent = ["CFD","STRUCTURE"]) # Use multiple capsIntent
# 
#         infoDict = myAnalysis.getAnalysisInfo(printinfo = False,  infoDict = True)
#         self.assertEqual(infoDict["intent"], "CFD;STRUCTURE")
# 
# #=============================================================================-
#     # Adding/getting attributes
#     def test_attributes(self):
# 
#         # Check list attributes
#         self.myAnalysis.addAttribute("testAttr", [1, 2, 3])
#         self.assertEqual(self.myAnalysis.getAttribute("testAttr"), [1,2,3])
# 
#         # Check float attributes
#         self.myAnalysis.addAttribute("testAttr_2", 10.0)
#         self.assertEqual(self.myAnalysis.getAttribute("testAttr_2"), 10.0)
# 
#         # Check string attributes
#         self.myAnalysis.addAttribute("testAttr_3", "anotherAttribute")
#         self.assertEqual(self.myAnalysis.getAttribute("testAttr_3"), "anotherAttribute")
# 
#         # Check over writing attribute
#         self.myAnalysis.addAttribute("testAttr_2", 30.0)
#         self.assertEqual(self.myAnalysis.getAttribute("testAttr_2"), 30.0)
# 
#         #self.myAnalysis.addAttribute("arrayOfStrings", ["1", "2", "3", "4"])
#         #self.assertEqual(self.myAnalysis.getAttribute("arrayOfStrings"), ["1", "2", "3", "4"])
#     
#=============================================================================-
    # Save geometry
    def test_writeGeometry(self):
         
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1
 
        analysis = problem.makeAnalysis("fun3dAIM", name = self.analysisDir + str(self.iDir)); self.__class__.iDir += 1

        analysis.writeGeometry(self.fileWrite)
        self.assertTrue(os.path.isfile(self.fileWrite))

#=============================================================================-
#         # Test bounding box
#     def test_getBoundingBox(self):
#         
#         box = self.myAnalysis.getBoundingBox()
#        
#         self.assertAlmostEqual(box["Farfield"][0], -80, 3)
#         self.assertAlmostEqual(box["Farfield"][3],  80, 3)

if __name__ == '__main__':
    unittest.main()
