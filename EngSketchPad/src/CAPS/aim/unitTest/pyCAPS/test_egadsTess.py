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

#==============================================================================
    def test_invalid_Mesh_Lenght_Scale(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        # Set project name so a mesh file is generated
        myAnalysis.input.Proj_Name = "pyCAPS_egadsTess_Test"
        
        # Run silent
        myAnalysis.input.Mesh_Quiet_Flag = True

        # Set output grid format since a project name is being supplied - Tecplot  file
        myAnalysis.input.Mesh_Format = "Tecplot"

        # Set output grid format since a project name is being supplied - Tecplot  file
        myAnalysis.input.Mesh_Length_Factor = -1

        with self.assertRaises(pyCAPS.CAPSError) as e:
            myAnalysis.runAnalysis()

        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    def test_setInput(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        myAnalysis.input.Proj_Name = "pyCAPS_EGADS_Test"
        myAnalysis.input.Mesh_Quiet_Flag = True
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

#==============================================================================
    def test_SingleBody_AnalysisOutVal(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        # Run silent
        myAnalysis.input.Mesh_Quiet_Flag = True

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

#==============================================================================
    def test_MultiBody(self):

        file = os.path.join("..","csmData","cfdMultiBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")
        
        # Run silent
        myAnalysis.input.Mesh_Quiet_Flag = True

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield" : {"tessParams" : [0.3*80, 0.2*80, 30]}}
        myAnalysis.input.Mesh_Sizing = Mesh_Sizing

        # Run
        myAnalysis.runAnalysis()

#==============================================================================
    def test_reenter(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

        # Run silent
        myAnalysis.input.Mesh_Quiet_Flag = True

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield": {"tessParams" : [0.3*80, 0.2*80, 30]}}
        myAnalysis.input.Mesh_Sizing = Mesh_Sizing

        myAnalysis.input.Mesh_Length_Factor = 1

        # Run 1st time
        myAnalysis.runAnalysis()
        
        NumberOfElement_1 = myAnalysis.output.NumberOfElement

        myAnalysis.input.Mesh_Length_Factor = 2

        # Run 2nd time coarser
        myAnalysis.runAnalysis()

        NumberOfElement_2 = myAnalysis.output.NumberOfElement

        self.assertGreater(NumberOfElement_1, NumberOfElement_2)


#==============================================================================
    def test_box(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Box",
                                                   capsIntent = ["box", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

#==============================================================================
    def test_cylinder(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Cylinder",
                                                   capsIntent = ["cylinder", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

#==============================================================================
    def test_cone(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Cone",
                                                   capsIntent = ["cone", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

#==============================================================================
    def test_torus(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Torus",
                                                   capsIntent = ["torus", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.view()

#==============================================================================
    def test_sphere(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Sphere",
                                                   capsIntent = ["sphere", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.view()

#==============================================================================
    def test_boxhole(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "BoxHole",
                                                   capsIntent = ["boxhole", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.view()

#==============================================================================
    def test_bullet(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Bullet",
                                                   capsIntent = ["bullet", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.view()

#==============================================================================
    def test_all(self):

        # Load egadsTess aim
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "All",
                                                   capsIntent = ["box", "cylinder", "cone", "torus", "sphere", "boxhole", "bullet", "farfield"]) 

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        egadsTess.runAnalysis()

        #egadsTess.view()

#==============================================================================
    def test_phase(self):

        file = os.path.join("..","csmData","cornerGeom.csm")
        
        probemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase0", capsFile=file, outLevel=0)

        egadsTess = myProblem.analysis.create(aim = "egadsTessAIM",
                                              name = "egadsTess",
                                              capsIntent = ["box", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"box": {"tessParams" : [0.3, 0.2, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing

        egadsTess.input.Mesh_Length_Factor = 1
        
        NumberOfNode_1    = egadsTess.output.NumberOfNode
        NumberOfElement_1 = egadsTess.output.NumberOfElement

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase1", prevPhaseName="Phase0", outLevel=0)

        egadsTess = myProblem.analysis["egadsTess"]
        
        # Check that the same outputs are still available
        self.assertTrue(egadsTess.output.Done)
        self.assertEqual(NumberOfNode_1   , egadsTess.output.NumberOfNode   )
        self.assertEqual(NumberOfElement_1, egadsTess.output.NumberOfElement)

        # Coarsen the mesh
        egadsTess.input.Mesh_Length_Factor = 2
        
        NumberOfNode_2    = egadsTess.output.NumberOfNode
        NumberOfElement_2 = egadsTess.output.NumberOfElement

        # Check that the counts have decreased
        self.assertGreater(NumberOfNode_1   , NumberOfNode_2   )
        self.assertGreater(NumberOfElement_1, NumberOfElement_2)

#==============================================================================
    def run_journal(self, myProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load egadsAIM
        if verbose: print(6*"-", line,"Load egadsTessAIM")
        egadsTess = myProblem.analysis.create(aim = "egadsTessAIM",
                                              capsIntent = ["box", "farfield"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run silent
        if verbose: print(6*"-", line,"Modify Mesh_Quiet_Flag")
        egadsTess.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Modify local mesh sizing parameters
        if verbose: print(6*"-", line,"Modify Tess_Params")
        Mesh_Sizing = {"box": {"tessParams" : [0.3, 0.2, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Length_Factor")
        egadsTess.input.Mesh_Length_Factor = 1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run 1st time
        if verbose: print(6*"-", line,"egadsTess runAnalysis")
        egadsTess.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"egadsTess NumberOfNode_1")
        NumberOfNode_1    = egadsTess.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"egadsTess NumberOfElement_1")
        NumberOfElement_1 = egadsTess.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())


        if verbose: print(6*"-", line,"Modify Mesh_Length_Factor")
        egadsTess.input.Mesh_Length_Factor = 2; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run 2nd time coarser
        if verbose: print(6*"-", line,"egadsTess runAnalysis")
        egadsTess.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"egadsTess NumberOfNode_2")
        NumberOfNode_2    = egadsTess.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"egadsTess NumberOfElement_2")
        NumberOfElement_2 = egadsTess.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Check that the counts have decreased
        self.assertGreater(NumberOfNode_1   , NumberOfNode_2   )
        self.assertGreater(NumberOfElement_1, NumberOfElement_2)

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def test_journal(self):

        capsFile = os.path.join("..","csmData","cornerGeom.csm")
        problemName = self.problemName+str(self.iProb)
        
        myProblem = pyCAPS.Problem(problemName, capsFile=capsFile, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_journal(myProblem, -1)
        
        myProblem.close()
        shutil.rmtree(problemName)
        
        #print(80*"=")
        #print(80*"=")
        # Create the problem to start journaling
        myProblem = pyCAPS.Problem(problemName, capsFile=capsFile, outLevel=0)
        myProblem.close()
        
        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(problemName, capsFile=capsFile, outLevel=0, useJournal=True)
            self.run_journal(myProblem, line_exit)
            myProblem.close()
            
        self.__class__.iProb += 1
        
if __name__ == '__main__':
    unittest.main()
