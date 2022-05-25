from __future__ import print_function
import unittest

import os
import glob
import shutil

import pyCAPS

class TestTETGEN(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_tetgenTest"
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

#==============================================================================
    def test_setInput(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        tetgen = myProblem.analysis.create(aim = "tetgenAIM")

        tetgen.input.Proj_Name = "pyCAPS_TetGen_Test"
        tetgen.input.Preserve_Surf_Mesh = True
        tetgen.input.Mesh_Verbose_Flag = False
        tetgen.input.Mesh_Quiet_Flag = False
        tetgen.input.Quality_Rad_Edge = 1.5
        tetgen.input.Quality_Angle = 0.0
        tetgen.input.Mesh_Format = "Tecplot"
        tetgen.input.Mesh_Gen_Input_String = ""
        tetgen.input.Mesh_Tolerance = 1e-16
        tetgen.input.Multiple_Mesh = False

        tetgen.input.Regions = {
              'A': { 'id': 10, 'seed': [0, 0,  1] },
              'B': { 'id': 20, 'seed': [0, 0, -1] }
            }

        tetgen.input.Holes = {
              'A': { 'seed': [0, 0,  1] },
              'B': { 'seed': [0, 0, -1] }
            }

#==============================================================================
    def test_reenter(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        egadsTess = myProblem.analysis.create(aim = "egadsTessAIM")

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield": {"tessParams" : [0.3*80, 0.2*80, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing

        egadsTess.input.Proj_Name = "test"
        egadsTess.input.Mesh_Format = "Tecplot"

        tetgen = myProblem.analysis.create(aim = "tetgenAIM",
                                           name = "tetgen")

        tetgen.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        tetgen.input.Mesh_Quiet_Flag = True
        tetgen.input.Proj_Name = "test"
        tetgen.input.Mesh_Format = "Tecplot"

        # Explicitly trigger mesh generation
        tetgen.runAnalysis()

        tetgen.input.Quality_Rad_Edge = 2.0

        # Explicitly trigger mesh generation again
        tetgen.runAnalysis()

#==============================================================================
    def test_phase(self):
        file = os.path.join("..","csmData","cfdSingleBody.csm")
        
        probemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase0", capsFile=file, outLevel=0)

        egadsTess = myProblem.analysis.create(aim = "egadsTessAIM",
                                              name = "egadsTess")

        egadsTess.input.Mesh_Quiet_Flag = True

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield": {"tessParams" : [0.3*80, 0.2*80, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing

        tetgen = myProblem.analysis.create(aim = "tetgenAIM",
                                           name = "tetgen")

        tetgen.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        tetgen.input.Mesh_Quiet_Flag = True
        tetgen.input.Proj_Name = "test"
        tetgen.input.Mesh_Format = "Tecplot"

        VolNumberOfNode_1    = tetgen.output.NumberOfNode
        VolNumberOfElement_1 = tetgen.output.NumberOfElement

        SurfNumberOfNode_1    = egadsTess.output.NumberOfNode
        SurfNumberOfElement_1 = egadsTess.output.NumberOfElement

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase1", prevPhaseName="Phase0", outLevel=0)

        egadsTess = myProblem.analysis["egadsTess"]
        tetgen    = myProblem.analysis["tetgen"]

        # Check that the same outputs are still available
        self.assertEqual(VolNumberOfNode_1   , tetgen.output.NumberOfNode   )
        self.assertEqual(VolNumberOfElement_1, tetgen.output.NumberOfElement)

        self.assertEqual(SurfNumberOfNode_1   , egadsTess.output.NumberOfNode   )
        self.assertEqual(SurfNumberOfElement_1, egadsTess.output.NumberOfElement)

        # Coarsen the mesh
        egadsTess.input.Mesh_Length_Factor = 2

        VolNumberOfNode_2    = tetgen.output.NumberOfNode
        VolNumberOfElement_2 = tetgen.output.NumberOfElement

        SurfNumberOfNode_2    = egadsTess.output.NumberOfNode
        SurfNumberOfElement_2 = egadsTess.output.NumberOfElement

        # Check that the counts have decreased
        self.assertGreater(VolNumberOfNode_1   , VolNumberOfNode_2   )
        self.assertGreater(VolNumberOfElement_1, VolNumberOfElement_2)

        self.assertGreater(SurfNumberOfNode_1   , SurfNumberOfNode_2   )
        self.assertGreater(SurfNumberOfElement_1, SurfNumberOfElement_2)


#==============================================================================
    def run_journal(self, myProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load egadsAIM
        if verbose: print(6*"-", line, line,"Load egadsAIM")
        egadsTess = myProblem.analysis.create(aim = "egadsTessAIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line, line,"egadsTess Mesh_Quiet_Flag")
        egadsTess.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Modify local mesh sizing parameters
        if verbose: print(6*"-", line, line,"egadsTess Mesh_Sizing")
        Mesh_Sizing = {"Farfield": {"tessParams" : [0.3*80, 0.2*80, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line, line,"egadsTess Mesh_Length_Factor")
        egadsTess.input.Mesh_Length_Factor = 1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Create the tetgen AIM
        if verbose: print(6*"-", line, line,"Load tetgenAIM")
        tetgen = myProblem.analysis.create(aim = "tetgenAIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line, line,"tetgen Mesh_Quiet_Flag")
        tetgen.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Link the surface mesh
        if verbose: print(6*"-", line, line,"tetgen Link Surface_Mesh")
        tetgen.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run egads explicitly
        #if verbose: print(6*"-", line, line,"Run EGADS")
        #egadsTess.runAnalysis(); line += 1
        #if line == line_exit: return line
        #if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run tetgen explicitly
        if verbose: print(6*"-", line, line,"Run TetGen")
        tetgen.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line, line,"tetgen VolNumberOfNode_1")
        VolNumberOfNode_1    = tetgen.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        
        if verbose: print(6*"-", line, line,"tetgen VolNumberOfNode_1")
        VolNumberOfElement_1 = tetgen.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Coarsen the mesh
        if verbose: print(6*"-", line, line,"Modify Quality_Rad_Edge")
        tetgen.input.Quality_Rad_Edge = 2.0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line, line,"tetgen VolNumberOfNode_2")
        VolNumberOfNode_2    = tetgen.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        
        if verbose: print(6*"-", line, line,"tetgen VolNumberOfElement_2")
        VolNumberOfElement_2 = tetgen.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        
        # Check that the counts have decreased
        self.assertGreater(VolNumberOfNode_1   , VolNumberOfNode_2   )
        self.assertGreater(VolNumberOfElement_1, VolNumberOfElement_2)

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def test_journal(self):

        capsFile = os.path.join("..","csmData","cfdSingleBody.csm")
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
