import unittest

import os
import glob
import shutil

import pyCAPS

from Mesh_Formats import Mesh_Formats
Mesh_Formats = Mesh_Formats.copy()

class TestInterference(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_interferenceTest"
        cls.iProb = 1
        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
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

        file = os.path.join("..","csmData","interference.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        interference = capsProblem.analysis.create(aim = "interferenceAIM")

        interference.input.Attr_Name = "_name"
        interference.input.OML = True
        interference.input.Tess_Params = [0.1, 0.01, 15.0]

#==============================================================================
    def test_OML(self):
        file = os.path.join("..","csmData","interference.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        interference = capsProblem.analysis.create(aim = "interferenceAIM")

        interference.input.Attr_Name = "_name"
        interference.input.OML = True
        interference.input.Tess_Params = [0.01, 0.001, 15.0]

        NamesTrue = ['Box', 'Sphere', 'Cylinder', 'OuterBox']
        Names = interference.output.Names
        #print(Names)

        DistancesTrue = [[0.0, -0.021330732618004396, 0.020000000000000018, 0.0],
                         [-0.021330732618004396, 0.0, 0.8031034039722986, -0.02143181466040841],
                         [0.020000000000000018, 0.8031034039722986, 0.0, 0.0],
                         [0.0, -0.02143181466040841, 0.0, 1.0]]
        Distances = interference.output.Distances
        #print(Distances)

        VolumesTrue = [1.0, 0.5883056710477834, 0.7235322985623731, 40.0]
        Volumes = interference.output.Volumes
        #print(Volumes)

        AreasTrue = [6.0, 3.3958811605473826, 4.462690758020493, 72.0]
        Areas = interference.output.Areas
        #print(Areas)

        CGsTrue = [[0.5, 0.5, 0.5],
                   [1.5, 0.5, 0.5],
                   [0.0, 1.5, 0.5],
                   [0.0, 0.0, 0.25]]
        CGs = interference.output.CGs
        #print(CGs)

        InertiasTrue = [[0.16666666666665342, -9.042245349977218e-06, -9.042245365631363e-06, -9.042245349977218e-06, 0.1666666666666422, -9.042245371237989e-06, -9.042245365631363e-06, -9.042245371237989e-06, 0.1666666666666149],
                        [0.06356676850723908, -4.2827996798999735e-07, -1.5283371085361708e-07, -4.2827996798999735e-07, 0.06356639663228969, 9.740635023847144e-07, -1.5283371085361708e-07, 9.740635023847144e-07, 0.06356707915949711],
                        [0.10194750753638093, 4.343775663947931e-08, -1.5418160573804762e-05, 4.343775663947931e-08, 0.101947507536407, 1.2212453270876722e-15, -1.5418160573804762e-05, 1.2212453270876722e-15, 0.08330629864571382],
                        [74.16666666656629, -7.233796296678933e-05, 4.262155949355529e-14, -7.233796296678933e-05, 74.16666666656629, -3.205676826420477e-14, 4.262155949355529e-14, -3.205676826420477e-14, 106.66666666675358]]
        Inertias = interference.output.Inertias
        #print(Inertias)

        for i in range(len(NamesTrue)):
            self.assertEqual(Names[i], NamesTrue[i])

        for i in range(len(DistancesTrue)):
            for j in range(len(DistancesTrue[i])):
                self.assertAlmostEqual(Distances[i][j], DistancesTrue[i][j], 2)

        for i in range(len(VolumesTrue)):
            self.assertAlmostEqual(Volumes[i], VolumesTrue[i], 3)

        for i in range(len(AreasTrue)):
            self.assertAlmostEqual(Areas[i], AreasTrue[i], 3)

        for i in range(len(CGsTrue)):
            for j in range(len(CGsTrue[i])):
                self.assertAlmostEqual(CGs[i][j], CGsTrue[i][j], 3)

        for i in range(len(InertiasTrue)):
            for j in range(len(InertiasTrue[i])):
                self.assertAlmostEqual(Inertias[i][j], InertiasTrue[i][j], 3)


#==============================================================================
    def test_boxes(self):
        file = os.path.join("..","csmData","interference_boxes.csm")
        capsProblem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=file, outLevel=0); self.__class__.iProb += 1

        interference = capsProblem.analysis.create(aim = "interferenceAIM")

        interference.input.Attr_Name = "_name"
        interference.input.OML = False
        interference.input.Tess_Params = [0.1, 0.001, 15.0]

        NamesTrue = ['Box1', 'Box2', 'Box3', 'Box4', 'Box5', 'Box6', 'Box7', 'Box8']
        Names = interference.output.Names
        #print(Names)

        DistancesTrue = [[0.0, 4.0, 1.0, 4.123105625617661, 4.0, 5.656854249492381, 7.0, 8.06225774829855],
                         [4.0, 0.0, 4.123105625617661, 1.0, 5.656854249492381, 4.0, 8.06225774829855, 7.0],
                         [1.0, 4.123105625617661, 0.0, 4.0, 1.0, 4.123105625617661, 4.0, 5.656854249492381],
                         [4.123105625617661, 1.0, 4.0, 0.0, 4.123105625617661, 1.0, 5.656854249492381, 4.0],
                         [4.0, 5.656854249492381, 1.0, 4.123105625617661, 0.0, 4.0, 1.0, 4.123105625617661],
                         [5.656854249492381, 4.0, 4.123105625617661, 1.0, 4.0, 0.0, 4.123105625617661, 1.0],
                         [7.0, 8.06225774829855, 4.0, 5.656854249492381, 1.0, 4.123105625617661, 0.0, 4.0],
                         [8.06225774829855, 7.0, 5.656854249492381, 4.0, 4.123105625617661, 1.0, 4.0, 0.0]]
        Distances = interference.output.Distances
        #print(Distances)

        VolumesTrue = [8.0]*8
        Volumes = interference.output.Volumes
        #print(Volumes)

        AreasTrue = [24.0]*8
        Areas = interference.output.Areas
        #print(Areas)

        CGsTrue = [[1.00, 1.00, 1.00],
                   [1.00, 7.00, 1.00],
                   [4.00, 1.00, 1.00],
                   [4.00, 7.00, 1.00],
                   [7.00, 1.00, 1.00],
                   [7.00, 7.00, 1.00],
                   [10.0, 1.00, 1.00],
                   [10.0, 7.00, 1.00]]
        CGs = interference.output.CGs
        #print(CGs)

        InertiasTrue = [[5.33333333333322, -0.00115740740722714, -0.0011574074073692486, -0.00115740740722714, 5.333333333333126, -0.0011574074074180984, -0.0011574074073692486, -0.0011574074074180984, 5.3333333333328845],
                        [5.33333333331592, -0.0011574074064029105, -0.008101851851845865, -0.0011574074064029105, 5.33333333333319, -0.0011574074074331975, -0.008101851851845865, -0.0011574074074331975, 5.3333333333252995],
                        [5.333333333333226, -0.0011574074067297602, -0.0011574074073052998, -0.0011574074067297602, 5.333333333334593, -0.004629629629640419, -0.0011574074073052998, -0.004629629629640419, 5.333333333334338],
                        [5.333333333317455, -0.0011574074044062854, -0.008101851851854747, -0.0011574074044062854, 5.333333333335105, -0.004629629629803844, -0.008101851851854747, -0.004629629629803844, 5.33333333332223],
                        [5.333333333333133, -0.0011574074062323803, -0.0011574074068789741, -0.0011574074062323803, 5.333333333312623, -0.00810185185179968, -0.0011574074068789741, -0.00810185185179968, 5.333333333312396],
                        [5.333333333312964, -0.0011574073972724364, -0.008101851851549213, -0.0011574073972724364, 5.33333333331359, -0.008101851851613162, -0.008101851851549213, -0.008101851851613162, 5.333333333318592],
                        [5.333333333333277, -0.0011574074060263229, -0.0011574074073621432, -0.0011574074060263229, 5.333333333343148, -0.01157407407410016, -0.0011574074073621432, -0.01157407407410016, 5.333333333342921],
                        [5.333333333319047, -0.00115740739965986, -0.008101851852131858, -0.00115740739965986, 5.3333333333445125, -0.011574074074488294, -0.008101851852131858, -0.011574074074488294, 5.333333333325527]]
        Inertias = interference.output.Inertias
        #print(Inertias)

        for i in range(len(NamesTrue)):
            self.assertEqual(Names[i], NamesTrue[i])

        for i in range(len(DistancesTrue)):
            for j in range(len(DistancesTrue[i])):
                self.assertAlmostEqual(Distances[i][j], DistancesTrue[i][j], 3)

        for i in range(len(VolumesTrue)):
            self.assertAlmostEqual(Volumes[i], VolumesTrue[i], 3)

        for i in range(len(AreasTrue)):
            self.assertAlmostEqual(Areas[i], AreasTrue[i], 3)

        for i in range(len(CGsTrue)):
            for j in range(len(CGsTrue[i])):
                self.assertAlmostEqual(CGs[i][j], CGsTrue[i][j], 3)

        for i in range(len(InertiasTrue)):
            for j in range(len(InertiasTrue[i])):
                self.assertAlmostEqual(Inertias[i][j], InertiasTrue[i][j], 3)


if __name__ == '__main__':
    unittest.main()
