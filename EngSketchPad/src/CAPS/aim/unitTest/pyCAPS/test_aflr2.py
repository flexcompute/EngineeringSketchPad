import unittest

import os
import glob
import shutil

import pyCAPS

from Mesh_Formats import Mesh_Formats
Mesh_Formats = Mesh_Formats.copy()

class TestAFLR2(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_aflr2Test"
        cls.iDir = 1
        cls.iProb = 1
        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove problem directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_setInput(self):

        file = os.path.join("..","csmData","cfd2D.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr2AIM")

        myAnalysis.input.Proj_Name = "pyCAPS_aflr2_Test"
        myAnalysis.input.Tess_Params = [.25,.01,15]
        myAnalysis.input.Mesh_Quiet_Flag = True
        myAnalysis.input.Mesh_Format = "Tecplot"
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

#==============================================================================
    def test_reenter(self):

        file = os.path.join("..","csmData","cfd2D.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        myAnalysis = myProblem.analysis.create(aim = "aflr2AIM")

        myAnalysis.input.Mesh_Quiet_Flag = True

        myAnalysis.input.Mesh_Format = Mesh_Formats

        airfoil = {"numEdgePoints" : 100, "edgeDistribution" : "Tanh", "initialNodeSpacing" : [0.001, 0.001]}
        myAnalysis.input.Mesh_Sizing = {"Airfoil"   : airfoil,
                                        "TunnelWall": {"numEdgePoints" : 50},
                                        "InFlow"    : {"numEdgePoints" : 25},
                                        "OutFlow"   : {"numEdgePoints" : 25},
                                        "2DSlice"   : {"tessParams" : [1.00, .01, 45]}}

        myAnalysis.input.Proj_Name = "pyCAPS_aflr2_Tri"

        # Run 1st time to make triangle mesh
        myAnalysis.runAnalysis()

        myAnalysis.input.Mesh_Gen_Input_String = "mquad=1 mpp=3"
        myAnalysis.input.Proj_Name = "pyCAPS_aflr2_Quad"

        # 'fast' file format does not support quad
        Mesh_Formats.remove("fast")
        myAnalysis.input.Mesh_Format = Mesh_Formats

        # Run 2nd time to make quad mesh
        myAnalysis.runAnalysis()


#==============================================================================
    def test_phase(self):

        file = os.path.join("..","csmData","cfd2D.csm")

        problemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile=file, outLevel=0)

        aflr2 = myProblem.analysis.create(aim = "aflr2AIM",
                                          name = "aflr2")

        # Run silent
        aflr2.input.Mesh_Quiet_Flag = True

        aflr2.input.Edge_Point_Min = 10
        aflr2.input.Edge_Point_Max = 30

        NumberOfNode_1    = aflr2.output.NumberOfNode
        NumberOfElement_1 = aflr2.output.NumberOfElement

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0)

        aflr2 = myProblem.analysis["aflr2"]

        # Check that the same outputs are still available
        self.assertEqual(NumberOfNode_1   , aflr2.output.NumberOfNode   )
        self.assertEqual(NumberOfElement_1, aflr2.output.NumberOfElement)

        # Coarsen the mesh
        aflr2.input.Edge_Point_Min = 5
        aflr2.input.Edge_Point_Max = 9

        NumberOfNode_2    = aflr2.output.NumberOfNode
        NumberOfElement_2 = aflr2.output.NumberOfElement

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
        if verbose: print(6*"-","Load aflr2AIM", line)
        aflr2 = myProblem.analysis.create(aim = "aflr2AIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run silent
        if verbose: print(6*"-","Modify Mesh_Quiet_Flag", line)
        aflr2.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Modify local mesh sizing parameters
        if verbose: print(6*"-","Modify Edge_Point_Min", line)
        aflr2.input.Edge_Point_Min = 10; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Modify Edge_Point_Min", line)
        aflr2.input.Edge_Point_Max = 30; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run 1st time
        if verbose: print(6*"-","aflr2 runAnalysis", line)
        aflr2.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr2 NumberOfNode_1", line)
        NumberOfNode_1    = aflr2.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr2 NumberOfElement_1", line)
        NumberOfElement_1 = aflr2.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run 2nd time coarser
        if verbose: print(6*"-","Modify Edge_Point_Min", line)
        aflr2.input.Edge_Point_Min = 5; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Modify Edge_Point_Min", line)
        aflr2.input.Edge_Point_Max = 20; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr2 NumberOfNode_2", line)
        NumberOfNode_2    = aflr2.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","aflr2 NumberOfElement_2", line)
        NumberOfElement_2 = aflr2.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Check that the counts have decreased
        self.assertGreater(NumberOfNode_1   , NumberOfNode_2   )
        self.assertGreater(NumberOfElement_1, NumberOfElement_2)

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def test_journal(self):

        capsFile = os.path.join("..","csmData","cfd2D.csm")
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
            myProblem = pyCAPS.Problem(problemName, phaseName="Scratch", capsFile=capsFile, outLevel=0)
            self.run_journal(myProblem, line_exit)
            myProblem.close()

        self.__class__.iProb += 1


if __name__ == '__main__':
    unittest.main()
