import unittest

import os
import glob
import shutil

import pyCAPS

from Mesh_Formats import Mesh_Formats
Mesh_Formats = Mesh_Formats.copy()

class TestCurveTess(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_curveTest"
        cls.iProb = 1
        cls.cleanUp()

        # Initialize a global Problem object
        cornerFile = os.path.join("..","csmData","cornerGeom.csm")
        cls.capsProblem = pyCAPS.Problem(cls.problemName, capsFile=cornerFile, outLevel=0)
        cls.egadsTess = cls.capsProblem.analysis.create(aim = "egadsTessAIM")

        # Run silent
        cls.egadsTess.input.Mesh_Quiet_Flag = True

    @classmethod
    def tearDownClass(cls):
        del cls.egadsTess
        del cls.capsProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_setInput(self):

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        
        curveTess.input["Surface_Mesh"].link(self.egadsTess.output["Surface_Mesh"])

        curveTess.input.Proj_Name = "pyCAPS_CurveTess_Test"
        curveTess.input.Element_Order = 2
        curveTess.input.Element_Class = 1
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input.Mesh_Format = "Exodus"
        curveTess.input.Mesh_Morph = False

#==============================================================================
    def test_invalid_inputs(self):

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")

        with self.assertRaises(pyCAPS.CAPSError) as e:
            curveTess.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")
        
        curveTess.input["Surface_Mesh"].link(self.egadsTess.output["Surface_Mesh"])

        curveTess.input.Element_Order = -1
        with self.assertRaises(pyCAPS.CAPSError) as e:
            curveTess.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        curveTess.input.Element_Order = 1
        curveTess.input.Element_Class = -1
        with self.assertRaises(pyCAPS.CAPSError) as e:
            curveTess.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        curveTess.input.Element_Order = 1
        curveTess.input.Element_Class = 3
        with self.assertRaises(pyCAPS.CAPSError) as e:
            curveTess.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        curveTess.input.Element_Order = 2
        curveTess.input.Element_Class = -1
        with self.assertRaises(pyCAPS.CAPSError) as e:
            curveTess.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        curveTess.input.Element_Order = 2
        curveTess.input.Element_Class = 3
        with self.assertRaises(pyCAPS.CAPSError) as e:
            curveTess.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    def test_reenter(self):

        file = os.path.join("..","csmData","cornerGeom.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        egadsTess = capsProblem.analysis.create(aim = "egadsTessAIM")
        egadsTess.input.Mesh_Quiet_Flag = True

        curveTess = capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Run 1st time
        egadsTess.input.Mesh_Length_Factor = 4

        NumberOfElement_1 = curveTess.output.NumberOfElement

        # Run 2nd time coarser
        egadsTess.input.Mesh_Length_Factor = 8

        NumberOfElement_2 = curveTess.output.NumberOfElement

        self.assertGreater(NumberOfElement_1, NumberOfElement_2)


#==============================================================================
    def test_box(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                     name = "Box",
                                                     capsIntent = ["box", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_cylinder(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                     name = "Cylinder",
                                                     capsIntent = ["cylinder", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_cone(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Cone",
                                                   capsIntent = ["cone", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_torus(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Torus",
                                                   capsIntent = ["torus", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_sphere(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Sphere",
                                                   capsIntent = ["sphere", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_boxhole(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "BoxHole",
                                                   capsIntent = ["boxhole", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_bullet(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "Bullet",
                                                   capsIntent = ["bullet", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_nodeBody(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "nodeBody",
                                                   capsIntent = ["nodeBody", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        curveTess.runAnalysis()
        #curveTess.geometry.view()

#==============================================================================
    def test_all(self):

        # Load egadsTess aim
        egadsTess = self.capsProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "All",
                                                   capsIntent = ["box", "cylinder", "cone", "torus", "sphere", "boxhole", "bullet", "nodeBody", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Mesh_Length_Factor = 8

        curveTess = self.capsProblem.analysis.create(aim = "curveTessAIM")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        # Just make sure it runs without errors...
        egadsTess.input.Multiple_Mesh = 'SingleDomain'
        curveTess.runAnalysis()
        egadsTess.input.Multiple_Mesh = 'MultiDomain'
        curveTess.runAnalysis()
        egadsTess.input.Multiple_Mesh = 'MultiFile'
        curveTess.runAnalysis()

        #curveTess.geometry.view()

#==============================================================================
    def test_phase(self):

        file = os.path.join("..","csmData","cornerGeom.csm")

        problemName = self.problemName + "_Phase"
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile=file, outLevel=0)

        egadsTess = capsProblem.analysis.create(aim = "egadsTessAIM",
                                                name = "egadsTess",
                                                capsIntent = ["box", "farfield"])

        # Run silent
        egadsTess.input.Mesh_Quiet_Flag = True

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"box": {"tessParams" : [0.3, 0.2, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing

        egadsTess.input.Mesh_Length_Factor = 4

        curveTess = capsProblem.analysis.create(aim = "curveTessAIM",
                                                name = "curveTess")
        curveTess.input.Mesh_Quiet_Flag = True
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        NumberOfNode_1    = curveTess.output.NumberOfNode
        NumberOfElement_1 = curveTess.output.NumberOfElement

        capsProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0)

        egadsTess = capsProblem.analysis["egadsTess"]
        curveTess = capsProblem.analysis["curveTess"]

        # Check that the same outputs are still available
        self.assertEqual(NumberOfNode_1   , curveTess.output.NumberOfNode   )
        self.assertEqual(NumberOfElement_1, curveTess.output.NumberOfElement)

        # Coarsen the mesh
        egadsTess.input.Mesh_Length_Factor = 8

        NumberOfNode_2    = curveTess.output.NumberOfNode
        NumberOfElement_2 = curveTess.output.NumberOfElement

        # Check that the counts have decreased
        self.assertGreater(NumberOfNode_1   , NumberOfNode_2   )
        self.assertGreater(NumberOfElement_1, NumberOfElement_2)

#==============================================================================
    def run_journal(self, capsProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Load egadsAIM
        if verbose: print(6*"-", line,"Load egadsTessAIM")
        egadsTess = capsProblem.analysis.create(aim = "egadsTessAIM",
                                                capsIntent = ["box", "farfield"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run silent
        if verbose: print(6*"-", line,"Modify Mesh_Quiet_Flag")
        egadsTess.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Modify local mesh sizing parameters
        if verbose: print(6*"-", line,"Modify Tess_Params")
        Mesh_Sizing = {"box": {"tessParams" : [0.3*8, 0.2*8, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Length_Factor")
        egadsTess.input.Mesh_Length_Factor = 4; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())


        # Load curveTessAIM
        if verbose: print(6*"-", line,"Load curveTessAIM")
        curveTess = capsProblem.analysis.create(aim = "curveTessAIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run silent
        if verbose: print(6*"-", line,"Modify Mesh_Quiet_Flag")
        curveTess.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Link surface meshes
        if verbose: print(6*"-", line,"Load egadsTessAIM")
        curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run 1st time
        if verbose: print(6*"-", line,"curveTess runAnalysis")
        curveTess.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"curveTess NumberOfNode_1")
        NumberOfNode_1    = curveTess.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"curveTess NumberOfElement_1")
        NumberOfElement_1 = curveTess.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())


        if verbose: print(6*"-", line,"Modify Mesh_Length_Factor")
        egadsTess.input.Mesh_Length_Factor = 8; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run 2nd time coarser
        if verbose: print(6*"-", line,"curveTess runAnalysis")
        curveTess.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"curveTess NumberOfNode_2")
        NumberOfNode_2    = curveTess.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"curveTess NumberOfElement_2")
        NumberOfElement_2 = curveTess.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Check that the counts have decreased
        self.assertGreater(NumberOfNode_1   , NumberOfNode_2   )
        self.assertGreater(NumberOfElement_1, NumberOfElement_2)

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def test_journal(self):

        capsFile = os.path.join("..","csmData","cornerGeom.csm")
        problemName = self.problemName+str(self.iProb)

        capsProblem = pyCAPS.Problem(problemName, capsFile=capsFile, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_journal(capsProblem, -1)

        capsProblem.close()
        shutil.rmtree(problemName)

        #print(80*"=")
        #print(80*"=")
        # Create the problem to start journaling
        capsProblem = pyCAPS.Problem(problemName, capsFile=capsFile, outLevel=0)
        capsProblem.close()

        for line_exit in range(line_total):
            #print(80*"=")
            capsProblem = pyCAPS.Problem(problemName, phaseName="Scratch", capsFile=capsFile, outLevel=0)
            self.run_journal(capsProblem, line_exit)
            capsProblem.close()

        self.__class__.iProb += 1

if __name__ == '__main__':
    unittest.main()
