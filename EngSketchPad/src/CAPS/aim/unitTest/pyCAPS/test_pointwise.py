import unittest
import os
import glob
import shutil

import platform
import time

# Import pyCAPS class file
import pyCAPS

# Helper function to check if an executable exists
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

class TestPointwise(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.problemName = "workDir_PointwiseTest"
        cls.iProb = 1

        cls.cleanUp()

        # Initialize Problem object
        cls.myProblem = pyCAPS.Problem(problemName = cls.problemName,
                                       capsFile = os.path.join("..","csmData","cornerGeom.csm"))

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

    # This is executed prior to each test
    def setUp(self):
        if platform.system() == "Windows":
            try:
                PW_INSTALL = os.environ["PW_HOME"]
            except:
                self.skipTest("No pointwise executable")
        else:
            if which("pointwise") == None:
                self.skipTest("No pointwise executable")

    def run_pointwise(self, pointwise):

        # Run AIM pre-analysis
        pointwise.preAnalysis()

        ####### Run pointwise ####################
        print ("\n\nRunning pointwise......")

        # Run pointwise via system call
        # try up to 30 times in case license is not available
        CAPS_GLYPH = os.environ["CAPS_GLYPH"]
        for i in range(30):
            try:
                if platform.system() == "Windows":
                    PW_HOME = os.environ["PW_HOME"]
                    pointwise.system(PW_HOME + r"\win64\bin\tclsh.exe " + CAPS_GLYPH + r"\GeomToMesh.glf caps.egads capsUserDefaults.glf")
                else:
                    pointwise.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
            except pyCAPS.CAPSError:
                time.sleep(10) # wait and try again
                continue

            time.sleep(1) # let the harddrive breathe
            if os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.gma')) and \
              (os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.ugrid')) or \
               os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.lb8.ugrid'))): break
            time.sleep(10) # wait and try again

        # make sure the execution was successful
        self.assertTrue( os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.gma')) and
                        (os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.ugrid')) or
                         os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.lb8.ugrid'))))

        # Run AIM post-analysis
        pointwise.postAnalysis()

    def test_executeError(self):

        # Load pointwise aim
        self.pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                        name = "BoxError",
                                                        capsIntent = ["box", "farfield"])
        # Run AIM post-analysis
        self.pointwise.preAnalysis()

        # Check for error when pointwise was not executed
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.pointwise.postAnalysis()

        self.assertEqual(e.exception.errorName, "CAPS_IOERR")
        del self.pointwise

    def test_SingleBody(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb),
                                   capsFile = file); self.__class__.iProb += 1

        pointwise = myProblem.analysis.create(aim = "pointwiseAIM",
                                              name = "SingleBody")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 10
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 12

        # Run
        self.run_pointwise(pointwise)

    def test_MultiBody(self):

        file = os.path.join("..","csmData","cfdMultiBody.csm")
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb),
                                   capsFile = file); self.__class__.iProb += 1

        pointwise = myProblem.analysis.create(aim = "pointwiseAIM",
                                              name = "MultiBody")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 10
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 12

        # Run
        self.run_pointwise(pointwise)

    def test_reenter(self):

        file = "../csmData/cfdSingleBody.csm"
        myProblem = pyCAPS.Problem(self.problemName+str(self.iProb),
                                   capsFile = file); self.__class__.iProb += 1

        pointwise = myProblem.analysis.create(aim = "pointwiseAIM",
                                              name = "SingleBody2")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 10
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 12

        # Run 1st time
        self.run_pointwise(pointwise)

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Run 2nd time coarser
        self.run_pointwise(pointwise)


    def test_box(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "Box",
                                                   capsIntent = ["box", "farfield"])

        #pointwise.saveGeometry("Box")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)

    def test_cylinder(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "Cylinder",
                                                   capsIntent = ["cylinder", "farfield"])

        #pointwise.saveGeometry("Cylinder")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)


    def test_cone(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "Cone",
                                                   capsIntent = ["cone", "farfield"])

        #pointwise.saveGeometry("Cone")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)

    def off_test_torus(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "Torus",
                                                   capsIntent = ["torus", "farfield"])

        #pointwise.saveGeometry("Torus")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)

    def test_sphere(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "Sphere",
                                                   capsIntent = ["sphere", "farfield"])

        #pointwise.saveGeometry("Sphere")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)

        #pointwise.view()

    def off_test_boxhole(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "BoxHole",
                                                   capsIntent = ["boxhole", "farfield"])

        #pointwise.saveGeometry("BoxHole")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)

    def test_bullet(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "Bullet",
                                                   capsIntent = ["bullet", "farfield"])

        #pointwise.saveGeometry("Bullet")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)

    def test_all(self):

        # Load pointwise aim
        pointwise = self.myProblem.analysis.create(aim = "pointwiseAIM",
                                                   name = "All",
                                                   capsIntent = ["box", "cone", "sphere", "bullet", "farfield"]) #, "cylinder", "torus", "boxhole"

        #pointwise.saveGeometry("GeomAll")

        # Global Min/Max number of points on edges
        pointwise.input.Connector_Initial_Dim = 6
        pointwise.input.Connector_Min_Dim     = 3
        pointwise.input.Connector_Max_Dim     = 8

        # Just make sure it runs without errors...
        self.run_pointwise(pointwise)

        #pointwise.view()

if __name__ == '__main__':
    unittest.main()
