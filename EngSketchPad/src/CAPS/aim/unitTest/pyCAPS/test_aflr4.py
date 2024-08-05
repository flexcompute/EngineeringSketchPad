import unittest

import os
import glob
import shutil

import pyCAPS

from Mesh_Formats import Mesh_Formats
Mesh_Formats = Mesh_Formats.copy()

class TestAFLR4(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_aflr4Test"
        cls.iProb = 1
        cls.cleanUp()

        # Initialize a global Problem object
        cornerFile = os.path.join("..","csmData","cornerGeom.csm")
        cls.capsProblem = pyCAPS.Problem(cls.problemName+str(cls.iProb), capsFile=cornerFile, outLevel=0); cls.iProb += 1

        # Single body
        file = os.path.join("..","csmData","cfdSingleBody.csm")
        cls.cfdSingleBody = pyCAPS.Problem(cls.problemName+str(cls.iProb), capsFile=file, outLevel=0); cls.iProb += 1

        cls.aflr4_sb = cls.cfdSingleBody.analysis.create(aim = "aflr4AIM")

        # Generate reference default mesh
        cls.aflr4_sb.input.Mesh_Quiet_Flag = True
        cls.aflr4_sb.runAnalysis()

    @classmethod
    def tearDownClass(cls):
        del cls.capsProblem
        del cls.aflr4_sb
        del cls.cfdSingleBody
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove problem directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_setInput(self):

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")

        aflr4.input.Proj_Name = "pyCAPS_AFLR4_Test"
        aflr4.input.Mesh_Quiet_Flag = True
        aflr4.input.Mesh_Format = "tecplot"
        aflr4.input.Mesh_Gen_Input_String = "auto_mode=0"
        aflr4.input.ff_cdfr = 1.4
        aflr4.input.min_ncell = 20
        aflr4.input.no_prox = True
        aflr4.input.abs_min_scale = 0.01
        aflr4.input.BL_thickness = 0.01
        aflr4.input.curv_factor = 0.4
        aflr4.input.max_scale = 0.1
        aflr4.input.min_scale = 0.01
        aflr4.input.Mesh_Length_Factor = 1.05
        aflr4.input.mer_all = True
        aflr4.input.erw_all = 0.7
        aflr4.input.Multiple_Mesh = "MultiDomain"
        aflr4.input.EGADS_Quad = False
        aflr4.input.AFLR4_Quad = False

#==============================================================================
    def test_invalid_Mesh_Lenght_Scale(self):

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")
        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = -1

        with self.assertRaises(pyCAPS.CAPSError) as e:
            aflr4.runAnalysis()

        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    def test_reenter(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 1

        aflr4.input.Mesh_Format = Mesh_Formats

        # Run 1st time
        aflr4.runAnalysis()

        numNodes = aflr4.output.NumberOfNode
        numElements = aflr4.output.NumberOfElement
        self.assertGreater(numNodes, 0)
        self.assertGreater(numElements, 0)

        aflr4.input.Mesh_Length_Factor = 2

        # Run 2nd time coarser
        aflr4.runAnalysis()

        self.assertLess(aflr4.output.NumberOfNode, numNodes)
        self.assertLess(aflr4.output.NumberOfElement, numElements)

#==============================================================================
    def test_Edge_scaleFactor(self):

        numNodes    = self.aflr4_sb.output.NumberOfNode
        numElements = self.aflr4_sb.output.NumberOfElement
        self.assertGreater(numNodes, 0)
        self.assertGreater(numElements, 0)

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")
        aflr4.input.Mesh_Quiet_Flag = True

        # Set mesh sizing parmeters
        aflr4.input.Mesh_Sizing = {"trailingEdge": {"scaleFactor":0.5}}

        self.assertGreater(aflr4.output.NumberOfNode, numNodes)
        self.assertGreater(aflr4.output.NumberOfElement, numElements)

#==============================================================================
    def test_Face_scaleFactor(self):

        numNodes    = self.aflr4_sb.output.NumberOfNode
        numElements = self.aflr4_sb.output.NumberOfElement
        self.assertGreater(numNodes, 0)
        self.assertGreater(numElements, 0)

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")
        aflr4.input.Mesh_Quiet_Flag = True

        # Set mesh sizing parmeters
        aflr4.input.Mesh_Sizing = {"leftWing": {"scaleFactor":0.5}}

        self.assertGreater(aflr4.output.NumberOfNode, numNodes)
        self.assertGreater(aflr4.output.NumberOfElement, numElements)

#==============================================================================
    def test_erw(self):

        numNodes    = self.aflr4_sb.output.NumberOfNode
        numElements = self.aflr4_sb.output.NumberOfElement
        self.assertGreater(numNodes, 0)
        self.assertGreater(numElements, 0)

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")
        aflr4.input.Mesh_Quiet_Flag = True

        # Increase mesh resolution only on left wing
        aflr4.input.Mesh_Sizing = {"leftWing": {"edgeWeight":1.0}}

        numNodes_erw    = aflr4.output.NumberOfNode
        numElements_erw = aflr4.output.NumberOfElement

        self.assertGreater(numNodes_erw   , numNodes)
        self.assertGreater(numElements_erw, numElements)

        aflr4.input.Mesh_Sizing = None

        # Increase mesh resolution by refining sharp edges
        aflr4.input.erw_all = 1.0
        aflr4.input.mer_all = True

        numNodes_erw_all    = aflr4.output.NumberOfNode
        numElements_erw_all = aflr4.output.NumberOfElement

        self.assertGreater(numNodes_erw_all   , numNodes_erw)
        self.assertGreater(numElements_erw_all, numElements_erw)

#==============================================================================
    def test_ff_cdfr(self):

        numNodes    = self.aflr4_sb.output.NumberOfNode
        numElements = self.aflr4_sb.output.NumberOfElement
        self.assertGreater(numNodes, 0)
        self.assertGreater(numElements, 0)

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")
        aflr4.input.Mesh_Quiet_Flag = True

        # Increase farfield mesh resolution
        aflr4.input.ff_cdfr = 1.1

        numNodes_ff    = aflr4.output.NumberOfNode
        numElements_ff = aflr4.output.NumberOfElement

        self.assertGreater(numNodes_ff   , numNodes)
        self.assertGreater(numElements_ff, numElements)

#==============================================================================
    def test_max_min_scale(self):

        numNodes    = self.aflr4_sb.output.NumberOfNode
        numElements = self.aflr4_sb.output.NumberOfElement
        self.assertGreater(numNodes, 0)
        self.assertGreater(numElements, 0)

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")
        aflr4.input.Mesh_Quiet_Flag = True

        # Coarsen mesh
        aflr4.input.max_scale = aflr4.input.max_scale*4
        aflr4.input.min_scale = aflr4.input.min_scale*4

        numNodes_scale    = aflr4.output.NumberOfNode
        numElements_scale = aflr4.output.NumberOfElement

        self.assertLess(numNodes_scale   , numNodes)
        self.assertLess(numElements_scale, numElements)

#==============================================================================
    def test_quad(self):

        numNodes    = self.aflr4_sb.output.NumberOfNode
        numElements = self.aflr4_sb.output.NumberOfElement
        numTri      = self.aflr4_sb.output.NumberOfTri
        numQuad     = self.aflr4_sb.output.NumberOfQuad
        self.assertGreater(numNodes   , 0)
        self.assertGreater(numElements, 0)
        self.assertGreater(numTri     , 0)
        self.assertEqual  (numQuad    , 0)

        aflr4 = self.cfdSingleBody.analysis.create(aim = "aflr4AIM")
        aflr4.input.Mesh_Quiet_Flag = True

        Quad_Mesh_Formats = Mesh_Formats.copy()
        if "fast"   in Quad_Mesh_Formats: Quad_Mesh_Formats.remove("fast")
        if "exodus" in Quad_Mesh_Formats: Quad_Mesh_Formats.remove("exodus")
        aflr4.input.Mesh_Format = Quad_Mesh_Formats

        aflr4.input.AFLR4_Quad = True

        # Run
        aflr4.runAnalysis()

        # Assert AnalysisOutVals
        self.assertNotEqual(aflr4.output.NumberOfNode   , numNodes)
        self.assertNotEqual(aflr4.output.NumberOfElement, numElements)
        self.assertNotEqual(aflr4.output.NumberOfTri    , numTri)
        self.assertGreater (aflr4.output.NumberOfQuad   , 0)

#==============================================================================
    def test_proximity(self):

        file = os.path.join("..","csmData","cfdMultiBody.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Sizing = {"Wing1": {"scaleFactor":2},
                                   "Wing2": {"scaleFactor":0.5}}

        # Run
        aflr4.runAnalysis()

        numNodes    = aflr4.output.NumberOfNode
        numElements = aflr4.output.NumberOfElement
        self.assertGreater(numNodes   , 0)
        self.assertGreater(numElements, 0)

        # Move the wings close to increase the resolution on wing1
        capsProblem.geometry.despmtr.winggap = 0.1

        # Run
        aflr4.runAnalysis()

        self.assertNotEqual(aflr4.output.NumberOfNode   , numNodes)
        self.assertNotEqual(aflr4.output.NumberOfElement, numElements)

#==============================================================================
    def test_regions(self):

        file = os.path.join("..","csmData","regions.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 1

        aflr4.input.Mesh_Format = Mesh_Formats

        # Run
        aflr4.runAnalysis()

#==============================================================================
    def test_skin(self):

        file = os.path.join("..","csmData","feaWingBEMAero.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True
        aflr4.input.skin = True

        aflr4.input.Mesh_Length_Factor = 1

        aflr4.input.Mesh_Format = Mesh_Formats

        # Run
        aflr4.runAnalysis()

#==============================================================================
    def test_box(self):

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
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
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
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
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
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
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
                                                 name = "Torus",
                                                 capsIntent = ["torus", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #self.capsProblem.geometry.save("torus.egads")
        # Just make sure it runs without errors...
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_sphere(self):

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
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
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
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
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
                                                 name = "Bullet",
                                                 capsIntent = ["bullet", "farfield"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #self.capsProblem.geometry.save("bullet.egads")
        # Just make sure it runs without errors...
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_all(self):

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM",
                                                 name = "All",
                                                 capsIntent = ["box", "cylinder", "cone", "torus", "sphere", "farfield", "bullet", "boxhole"])

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_faceMatch(self):

        file = os.path.join("..","csmData","multiRegions.csm")

        problemName = self.problemName + "_faceMatch"
        capsProblem = pyCAPS.Problem(problemName, capsFile=file, outLevel=0)

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4_2",
                                            capsIntent = ["regions2"])
        aflr4.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4_3",
                                            capsIntent = ["regions2", "regions3"])
        aflr4.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4_4",
                                            capsIntent = ["regions2", "regions3", "regions4"])
        aflr4.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4_5",
                                            capsIntent = ["regions2", "regions3", "regions4", "regions5"])
        aflr4.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4_6",
                                            capsIntent = ["regions2", "regions3", "regions4", "regions5", "regions6"])
        aflr4.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4_7",
                                            capsIntent = ["regions2", "regions3", "regions4", "regions5", "regions6", "regions7"])
        aflr4.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4_8",
                                            capsIntent = ["regions2", "regions3", "regions4", "regions5", "regions6", "regions7", "regions8"])
        aflr4.input.Mesh_Quiet_Flag = True

        # Just make sure it runs without errors...
        aflr4.input.Multiple_Mesh = 'SingleDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiDomain'
        aflr4.runAnalysis()
        aflr4.input.Multiple_Mesh = 'MultiFile'
        aflr4.runAnalysis()

        #aflr4.geometry.view()

#==============================================================================
    def test_faceMatch_MultiDomain(self):

        for csm in [
                    "example1",
                    "example2",
                    "example3",
                    "example4",
                   # "example5",
                    "example6",
                    "example7",
                    "example8",
                    "example9",
                    "multi_prim",
                    "cyl_seam",
                    ]:
            file = os.path.join("..","csmData","MultiDomain",csm+".csm")

            print(file)
            problemName = self.problemName + "_faceMatch_" + csm
            capsProblem = pyCAPS.Problem(problemName, capsFile=file, outLevel=0)

            aflr4 = capsProblem.analysis.create(aim = "aflr4AIM")

            aflr4.input.Mesh_Format = "tecplot"
            aflr4.input.Mesh_Quiet_Flag = True
            aflr4.input.Mesh_Length_Factor = 1
            aflr4.input.min_scale = 0.05
            #aflr4.input.max_scale = 0.5
            aflr4.input.Multiple_Mesh = 'SingleDomain'
            aflr4.runAnalysis()
            aflr4.input.Multiple_Mesh = 'MultiDomain'
            aflr4.runAnalysis()
            aflr4.input.Multiple_Mesh = 'MultiFile'
            aflr4.runAnalysis()

#==============================================================================
    def test_phase(self):

        file = os.path.join("..","csmData","cornerGeom.csm")

        problemName = self.problemName + "_Phase"
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile=file, outLevel=0)

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4",
                                            capsIntent = ["box", "farfield"])

        # Run silent
        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 1

        NumberOfNode_1    = aflr4.output.NumberOfNode
        NumberOfElement_1 = aflr4.output.NumberOfElement

        capsProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0)

        aflr4 = capsProblem.analysis["aflr4"]

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
    def run_journal(self, capsProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Load egadsAIM
        if verbose: print(6*"-","Load aflr4AIM", line)
        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            capsIntent = ["box", "farfield"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run silent
        if verbose: print(6*"-","Modify Mesh_Quiet_Flag", line)
        aflr4.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Modify local mesh sizing parameters
        if verbose: print(6*"-","Modify Mesh_Length_Factor", line)
        aflr4.input.Mesh_Length_Factor = 1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run 1st time
        if verbose: print(6*"-","aflr4 runAnalysis", line)
        aflr4.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfNode_1", line)
        NumberOfNode_1    = aflr4.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfElement_1", line)
        NumberOfElement_1 = aflr4.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run 2nd time coarser
        if verbose: print(6*"-","Modify Mesh_Length_Factor", line)
        aflr4.input.Mesh_Length_Factor = 2; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfNode_2", line)
        NumberOfNode_2    = aflr4.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-","aflr4 NumberOfElement_2", line)
        NumberOfElement_2 = aflr4.output.NumberOfElement; line += 1
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
