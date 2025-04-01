import unittest

import os
import glob
import shutil

import pyCAPS

from Mesh_Formats import Mesh_Formats
Mesh_Formats = Mesh_Formats.copy()

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

class TestMETRIS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_metrisTest"
        cls.iProb = 1
        cls.cleanUp()

        # Initialize a global Problem object
        cornerFile = os.path.join("..","csmData","cornerGeom.csm")
        cls.capsProblem = pyCAPS.Problem(cls.problemName, capsFile=cornerFile, outLevel=0)

    @classmethod
    def tearDownClass(cls):
        del cls.capsProblem
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

    # This is executed prior to each test
    def setUp(self):
        if which("ref") == None:
            self.skipTest("No 'ref' executable")
        if which("metris") == None:
            self.skipTest("No 'metris' executable")

#==============================================================================
    # Test inputs
    def test_inputs(self):

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input.metris = "metris"
        metris.input.Passes = 1
        metris.input.Mesh_Format = "Tecplot"
        metris.input.MetricFieldFile = os.path.join(metris.analysisDir, "metric.sol")

#==============================================================================
    def test_Area_Mesh(self):

        # Load CSM file
        capsProblem = pyCAPS.Problem(problemName = self.problemName + "2D",
                                     capsFile = os.path.join("..","csmData","cfd2D.csm"),
                                     outLevel = 0)

        # Put geometry in x-y plane
        capsProblem.geometry.cfgpmtr.plane = 2

        # Load aflr2 aim
        aflr2 = capsProblem.analysis.create(aim = "aflr2AIM")

        # Run silent
        aflr2.input.Mesh_Quiet_Flag = True

        aflr2.input.Edge_Point_Min = 5
        aflr2.input.Edge_Point_Max = 10

        # Load metris aim
        metris = capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr2.output["Area_Mesh"])

        # 'stl' file format does not support 2D
        Mesh_Formats.remove("stl")
        Mesh_Formats.remove("bstl")
        metris.input.Mesh_Format = Mesh_Formats

        # Load refine aim for MultiScale metric
        refine = capsProblem.analysis.create(aim = "refineAIM")

        refine.input["Mesh"].link(metris.output["Mesh"])

        xyz = refine.output.xyz

        refine.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = (xyz[i][0]**2 + xyz[i][1]**2)**(1/2.)

        with open(refine.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("2\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        refine.input.Complexity = 700

        metris.input.MetricFieldFile = refine.output.MetricFieldFile
        print("Run metris 1")
        metris.input["Mesh"].unlink()
        metris.runAnalysis()

        xyz = metris.output.xyz
    
        # Make sure the mesh changed
        self.assertNotEqual(len(xyz), len(scalar))

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = (xyz[i][0]**2 + xyz[i][1]**2)**(1/2.)

        with open(refine.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("2\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        refine.input.Complexity = refine.input.Complexity*1.5
        refine.input.Passes = 1
 
        metris.input.MetricFieldFile = refine.output.MetricFieldFile
        print("Run metris 2")
        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

        # refine.input.ScalarFieldFile = None
        # refine.input.HessianFieldFile = os.path.join(metris.analysisDir, "hessian.sol")
        #
        # hessian = [[0,0,0]]*len(xyz)
        # for i in range(len(xyz)):
        #     x = xyz[i][0]
        #     y = xyz[i][1]
        #     r = (x**2 + y**2)**(1/2.)
        #     r = max(r,1e-10)
        #     hessian[i][0] = 1/r - x**2/r**3 # d2r/dx2
        #     hessian[i][1] =     - x*y /r**3 # d2r/dxy
        #     hessian[i][2] = 1/r - y**2/r**3 # d2r/dy2
        #
        # with open(refine.input.HessianFieldFile, "w") as f:
        #     f.write("MeshVersionFormatted 2\n")
        #     f.write("Dimension\n")
        #     f.write("2\n")
        #     f.write("SolAtVertices\n")
        #     f.write(str(len(hessian)) + " 1 3\n")
        #     for i in range(len(hessian)):
        #         for j in range(3):
        #             f.write(str(hessian[i][j]) + " ")
        #         f.write("\n")
        #
        # refine.input.Complexity = 700
        # refine.input.Passes = 1
        #
        # metris.input.MetricFieldFile = refine.output.MetricFieldFile
        #
        # print("Run metris 3")
        # metris.runAnalysis()
        #
        # xyz = metris.output.xyz
        #
        # # Make sure the mesh changed
        # self.assertNotEqual(len(xyz), len(hessian))

#==============================================================================
    def off_test_cfdSingleBody(self):

        file = os.path.join("..","csmData","cfdSingleBody.csm")

        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                            name = "aflr4")

        # Modify local mesh sizing parameters
        aflr4.input.Mesh_Quiet_Flag = True

        #aflr4.input.max_scale = 2
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4
        aflr4.input.curv_factor = 2.0
        aflr4.input.erw_all = 0.0
        aflr4.input.mer_all = False

        aflr4.input.Mesh_Length_Factor = 4

        aflr3 = capsProblem.analysis.create(aim = "aflr3AIM",
                                            name = "aflr3")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set project name
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = capsProblem.analysis.create(aim = "metrisAIM",
                                             name = "metris")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])
        metris.input.Passes = 1 #Just for testing

        metris.input.Mesh_Format = Mesh_Formats

        su2 = capsProblem.analysis.create(aim = "su2AIM")

        su2.input["Mesh"].link(metris.output["Mesh"])

        su2.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"},
                                        "Farfield":"farfield"}

        scalar_last = []
        for iadapt in range(3):
            su2.preAnalysis()
            su2.postAnalysis()

            xyz = metris.output.xyz

            metris.input["Mesh"].unlink()
            metris.input.ScalarFieldFile = os.path.join(su2.analysisDir, "scalar.sol")

            scalar = [1]*len(xyz)

            with open(metris.input.ScalarFieldFile, "w") as f:
                f.write("Dimension\n")
                f.write("3\n")
                f.write("SolAtVertices\n")
                f.write(str(len(scalar)) + " 1 1\n")
                for i in range(len(xyz)):
                    f.write(str(scalar[i]) + "\n")

            metris.input.Complexity = len(scalar)/2

            self.assertNotEqual(len(scalar_last), len(scalar))
            scalar_last = scalar

#==============================================================================
    def off_test_box_volume(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 1
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 0
        self.capsProblem.geometry.cfgpmtr.boxhole  = 0
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/4
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

        metris.input.ScalarFieldFile = None
        metris.input.HessianFieldFile = os.path.join(metris.analysisDir, "hessian.sol")

        hessian = [[1,0,0,1,0,1]]*len(xyz)

        with open(metris.input.HessianFieldFile, "w") as f:
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(hessian)) + " 1 3\n")
            for i in range(len(hessian)):
                for j in range(6):
                    f.write(str(hessian[i][j]) + " ")
                f.write("\n")

        metris.input.Complexity = len(hessian)/4
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_fails_test_box_surface(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 1
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 0
        self.capsProblem.geometry.cfgpmtr.boxhole  = 0
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = False

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr4.output["Surface_Mesh"])

        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [1]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_fails_test_cylinder(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 0
        self.capsProblem.geometry.cfgpmtr.cylinder = 1
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 0
        self.capsProblem.geometry.cfgpmtr.boxhole  = 0
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_fails_test_cone(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 0
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 1
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 0
        self.capsProblem.geometry.cfgpmtr.boxhole  = 0
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_fails_test_torus(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 0
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 1
        self.capsProblem.geometry.cfgpmtr.sphere   = 0
        self.capsProblem.geometry.cfgpmtr.boxhole  = 0
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_test_sphere(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 0
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 1
        self.capsProblem.geometry.cfgpmtr.boxhole  = 0
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 2
        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_test_boxhole(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 0
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 0
        self.capsProblem.geometry.cfgpmtr.boxhole  = 1
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 2
        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_fails_test_bullet(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 0
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 0
        self.capsProblem.geometry.cfgpmtr.boxhole  = 0
        self.capsProblem.geometry.cfgpmtr.bullet   = 1
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 2
        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))

#==============================================================================
    def off_test_all(self):

        self.capsProblem.geometry.cfgpmtr.single   = 1
        self.capsProblem.geometry.cfgpmtr.box      = 1
        self.capsProblem.geometry.cfgpmtr.cylinder = 0
        self.capsProblem.geometry.cfgpmtr.cone     = 0
        self.capsProblem.geometry.cfgpmtr.torus    = 0
        self.capsProblem.geometry.cfgpmtr.sphere   = 1
        self.capsProblem.geometry.cfgpmtr.boxhole  = 1
        self.capsProblem.geometry.cfgpmtr.bullet   = 0
        self.capsProblem.geometry.cfgpmtr.nodebody = 0

        # Load aflr4 aim
        aflr4 = self.capsProblem.analysis.create(aim = "aflr4AIM")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 2
        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        #aflr4.geometry.save("box.egads")

        # Load aflr3 aim
        aflr3 = self.capsProblem.analysis.create(aim = "aflr3AIM")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        # Set output
        aflr3.input.Mesh_Quiet_Flag = True

        # Load metris aim
        metris = self.capsProblem.analysis.create(aim = "metrisAIM")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        metris.runAnalysis()

        xyz = metris.output.xyz

        self.assertNotEqual(len(xyz), len(scalar))


#==============================================================================
    def off_test_phase(self):
        file = os.path.join("..","csmData","cornerGeom.csm")

        problemName = self.problemName + "_Phase"
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile=file, outLevel=0)

        capsProblem.geometry.cfgpmtr.single   = 1
        capsProblem.geometry.cfgpmtr.box      = 1
        capsProblem.geometry.cfgpmtr.cylinder = 0
        capsProblem.geometry.cfgpmtr.cone     = 0
        capsProblem.geometry.cfgpmtr.torus    = 0
        capsProblem.geometry.cfgpmtr.sphere   = 0
        capsProblem.geometry.cfgpmtr.boxhole  = 0
        capsProblem.geometry.cfgpmtr.bullet   = 0
        capsProblem.geometry.cfgpmtr.nodebody = 0

        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM", name="aflr4")

        aflr4.input.Mesh_Quiet_Flag = True

        aflr4.input.Mesh_Length_Factor = 2
        aflr4.input.max_scale = 0.5
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr   = 1.4

        aflr3 = capsProblem.analysis.create(aim = "aflr3AIM", name="aflr3")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        aflr3.input.Mesh_Quiet_Flag = True

        VolNumberOfNode_1    = aflr3.output.NumberOfNode
        VolNumberOfElement_1 = aflr3.output.NumberOfElement

        SurfNumberOfNode_1    = aflr4.output.NumberOfNode
        SurfNumberOfElement_1 = aflr4.output.NumberOfElement

        # Load metris aim
        metris = capsProblem.analysis.create(aim = "metrisAIM", name="metris")

        metris.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        metris.input.Complexity = 50
        metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        xyz = metris.output.xyz

        metris.input["Mesh"].unlink()

        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2
        metris.input.Passes = 1

        xyz = metris.output.xyz

        RefNumberOfNode_1 = len(xyz)

        capsProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0)

        aflr4 = capsProblem.analysis["aflr4"]
        aflr3 = capsProblem.analysis["aflr3"]
        metris = capsProblem.analysis["metris"]

        # Check that the same outputs are still available
        self.assertEqual(VolNumberOfNode_1   , aflr3.output.NumberOfNode   )
        self.assertEqual(VolNumberOfElement_1, aflr3.output.NumberOfElement)

        self.assertEqual(SurfNumberOfNode_1   , aflr4.output.NumberOfNode   )
        self.assertEqual(SurfNumberOfElement_1, aflr4.output.NumberOfElement)

        xyz = metris.output.xyz

        self.assertEqual(RefNumberOfNode_1, len(xyz))

        metris.input.Passes = 1
        #metris.input.ScalarFieldFile = os.path.join(metris.analysisDir, "scalar.sol")

        # Coarsen the mesh
        scalar = [0]*len(xyz)
        for i in range(len(xyz)):
            scalar[i] = xyz[i][0]**2

        with open(metris.input.ScalarFieldFile, "w") as f:
            f.write("Dimension\n")
            f.write("3\n")
            f.write("SolAtVertices\n")
            f.write(str(len(scalar)) + " 1 1\n")
            for i in range(len(xyz)):
                f.write(str(scalar[i]) + "\n")

        metris.input.Complexity = len(scalar)/2

        xyz = metris.output.xyz

        RefNumberOfNode_2 = len(xyz)

        # Check that the counts have decreased
        self.assertNotEqual(RefNumberOfNode_1   , RefNumberOfNode_2   )

#==============================================================================
    def run_journal(self, capsProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Load egadsAIM
        if verbose: print(6*"-", line,"Load aflr4AIM")
        aflr4 = capsProblem.analysis.create(aim = "aflr4AIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Quiet_Flag")
        aflr4.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Length_Factor")
        aflr4.input.Mesh_Length_Factor = 20; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Create the aflr3 AIM
        if verbose: print(6*"-", line,"Load aflr3AIM")
        aflr3 = capsProblem.analysis.create(aim = "aflr3AIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Quiet_Flag")
        aflr3.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Link the surface mesh
        if verbose: print(6*"-", line,"Link Surface_Mesh")
        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run aflr4 explicitly
        #if verbose: print(6*"-", line,"Run AFLR4")
        #aflr4.runAnalysis(); line += 1
        #if line == line_exit: return line
        #if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Run aflr3 explicitly
        if verbose: print(6*"-", line,"Run AFLR3")
        aflr3.runAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"aflr3 VolNumberOfNode_1")
        VolNumberOfNode_1    = aflr3.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"aflr3 VolNumberOfNode_1")
        VolNumberOfElement_1 = aflr3.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())


        # Coarsen the mesh
        if verbose: print(6*"-", line,"Modify Mesh_Length_Factor")
        aflr4.input.Mesh_Length_Factor = 40; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"aflr3 VolNumberOfNode_2")
        VolNumberOfNode_2    = aflr3.output.NumberOfNode; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        if verbose: print(6*"-", line,"aflr3 VolNumberOfElement_2")
        VolNumberOfElement_2 = aflr3.output.NumberOfElement; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(capsProblem.journaling())

        # Check that the counts have decreased
        self.assertGreater(VolNumberOfNode_1   , VolNumberOfNode_2   )
        self.assertGreater(VolNumberOfElement_1, VolNumberOfElement_2)

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def off_test_journal(self):

        capsFile = os.path.join("..","csmData","cfdSingleBody.csm")
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
