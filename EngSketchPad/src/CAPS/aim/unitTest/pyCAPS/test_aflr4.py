from __future__ import print_function
import unittest

import os
import glob
import shutil

import pyCAPS

class TestAFLR4(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_aflr4Test"
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

        # Remove problem directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_invalid_Mesh_Lenght_Scale(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

        myAnalysis.input.Mesh_Quiet_Flag = True

        # Set project name so a mesh file is generated
        myAnalysis.input.Proj_Name = "pyCAPS_AFLR4_Test"

        # Set output grid format since a project name is being supplied - Tecplot  file
        myAnalysis.input.Mesh_Format = "Tecplot"

        myAnalysis.input.Mesh_Length_Factor = -1

        with self.assertRaises(pyCAPS.CAPSError) as e:
            myAnalysis.runAnalysis()

        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    def test_setInput(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

        myAnalysis.input.Proj_Name = "pyCAPS_AFLR4_Test"
        myAnalysis.input.Mesh_Quiet_Flag = True
        myAnalysis.input.Mesh_Format = "Tecplot"
        myAnalysis.input.Mesh_ASCII_Flag = True
        myAnalysis.input.Mesh_Gen_Input_String = "auto_mode=0"
        myAnalysis.input.ff_cdfr = 1.4
        myAnalysis.input.min_ncell = 20
        myAnalysis.input.no_prox = True
        myAnalysis.input.abs_min_scale = 0.01
        myAnalysis.input.BL_thickness = 0.01
        myAnalysis.input.curv_factor = 0.4
        myAnalysis.input.max_scale = 0.1
        myAnalysis.input.min_scale = 0.01
        myAnalysis.input.Mesh_Length_Factor = 1.05
        myAnalysis.input.erw_all = 0.7

#==============================================================================
    def test_SingleBody_AnalysisOutVal(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

        myAnalysis.input.Mesh_Quiet_Flag = True

        # Run
        myAnalysis.runAnalysis()

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

        myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

        myAnalysis.input.Mesh_Quiet_Flag = True

        # Run
        myAnalysis.runAnalysis()

#==============================================================================
    def test_reenter(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

        myAnalysis.input.Mesh_Quiet_Flag = True

        myAnalysis.input.Mesh_Length_Factor = 1

        # Run 1st time
        myAnalysis.runAnalysis()

        myAnalysis.input.Mesh_Length_Factor = 2

        # Run 2nd time coarser
        myAnalysis.runAnalysis()


#==============================================================================
    def test_box(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "Box",
                                               capsIntent = ["box", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Just make sure it runs without errors...
        aflr4.runAnalysis()

#==============================================================================
    def test_cylinder(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "Cylinder",
                                               capsIntent = ["cylinder", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Just make sure it runs without errors...
        aflr4.runAnalysis()

#==============================================================================
    def test_cone(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "Cone",
                                               capsIntent = ["cone", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Just make sure it runs without errors...
        aflr4.runAnalysis()

#==============================================================================
    def test_torus(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "Torus",
                                               capsIntent = ["torus", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #self.myProblem.geometry.save("torus.egads")
        # Just make sure it runs without errors...
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_sphere(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "Sphere",
                                               capsIntent = ["sphere", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Just make sure it runs without errors...
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_boxhole(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "BoxHole",
                                               capsIntent = ["boxhole", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Just make sure it runs without errors...
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_bullet(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "Bullet",
                                               capsIntent = ["bullet", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #self.myProblem.geometry.save("bullet.egads")
        # Just make sure it runs without errors...
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_all(self):

        # Load aflr4 aim
        aflr4 = self.myProblem.analysis.create(aim = "aflr4AIM",
                                               name = "All",
                                               capsIntent = ["box", "cylinder", "cone", "torus", "sphere", "farfield", "bullet", "boxhole"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Just make sure it runs without errors...
        aflr4.runAnalysis()

        #aflr4.geometry.view()


#==============================================================================
    def test_phase(self):

        file = os.path.join("..","csmData","cornerGeom.csm")
        
        probemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase0", capsFile=file, outLevel=0)

        aflr4 = myProblem.analysis.create(aim = "aflr4AIM",
                                          name = "aflr4",
                                          capsIntent = ["box", "farfield"])

        # Run silent
        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 1
        
        NumberOfNode_1    = aflr4.output.NumberOfNode
        NumberOfElement_1 = aflr4.output.NumberOfElement

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase1", prevPhaseName="Phase0", outLevel=0)

        aflr4 = myProblem.analysis["aflr4"]
        
        # Check that the same outputs are still available
        self.assertEqual(NumberOfNode_1   , aflr4.output.NumberOfNode   )
        self.assertEqual(NumberOfElement_1, aflr4.output.NumberOfElement)

        # Coarsen the mesh
        aflr4.input.Mesh_Length_Factor = 2
        
        NumberOfNode_2    = aflr4.output.NumberOfNode
        NumberOfElement_2 = aflr4.output.NumberOfElement

        # Check that the counts have decreased
        self.assertGreater(NumberOfNode_1   , NumberOfNode_2   )
        self.assertGreater(NumberOfElement_1, NumberOfElement_2)


#==============================================================================
    def run_journal(self, myProblem, line_exit):

        verbose = True

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load egadsAIM
        if verbose: print(6*"-","Load aflr4AIM", line)
        aflr4 = myProblem.analysis.create(aim = "aflr4AIM",
                                          capsIntent = ["box", "farfield"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run silent
        if verbose: print(6*"-","Modify Mesh_Quiet_Flag", line)
        aflr4.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Modify local mesh sizing parameters
        if verbose: print(6*"-","Modify Mesh_Length_Factor", line)
        aflr4.input.Mesh_Length_Factor = 1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run 1st time
        if verbose: print(6*"-","aflr4 runAnalysis", line)
        aflr4.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfNode_1", line)
        NumberOfNode_1    = aflr4.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfElement_1", line)
        NumberOfElement_1 = aflr4.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run 2nd time coarser
        if verbose: print(6*"-","Modify Mesh_Length_Factor", line)
        aflr4.input.Mesh_Length_Factor = 2; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfNode_2", line)
        NumberOfNode_2    = aflr4.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfElement_2", line)
        NumberOfElement_2 = aflr4.output.NumberOfElement; line += 1
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
