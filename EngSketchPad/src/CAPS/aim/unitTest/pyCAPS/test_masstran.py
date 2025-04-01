# Import other need modules
import unittest

import os
import glob
import shutil

# Import pyCAPS class file
import pyCAPS

class TestMasstran(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.problemName = "workDir_masstranTest"
        cls.iProb = 1
        cls.cleanUp()

        # Initialize Problem object
        file = os.path.join("..","csmData","masstran.csm")
        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=file, outLevel=0)

        # Load masstran aim
        cls.masstran = cls.myProblem.analysis.create(aim = "masstranAIM")

        # Min/Max number of points on an edge for quadding
        cls.masstran.input.Edge_Point_Min = 21
        cls.masstran.input.Edge_Point_Max = 21

        #cls.masstran.input.Tess_Params = [1/64, 0.001, 15.0]

        # Generate quad meshes
        cls.masstran.input.Quad_Mesh = True

        # Set materials
        madeupium    = {"materialType" : "isotropic",
                        "density"      : 10}
        unobtainium  = {"materialType" : "isotropic",
                        "density"      : 20}

        cls.masstran.input.Material = {"madeupium" : madeupium, "unobtainium" : unobtainium}

        # Set properties
        cls.shell1 = {"propertyType"      : "Shell",
                      "membraneThickness" : 2.0,
                      "material"          : "madeupium"}

        cls.shell2 = {"propertyType"      : "Shell",
                      "membraneThickness" : 3.0,
                      "material"          : "unobtainium"}

        cls.pointMass = {"propertyType" : "ConcentratedMass",
                         "mass"         : 3.0,
                         "massInertia"  : [10.0, 11.0, 12.0, 13.0, 14.0, 15.0]}

    @classmethod
    def tearDownClass(cls):
        del cls.myProblem
        del cls.masstran
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):
        # Remove problem directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def check_derivs(self, Itol = 5):

        self.masstran.input.Design_Variable = {"L":{}, "W":{}, "H":{}, 
                                               "x0":{"initialValue":0.0}, 
                                               "x1":{"initialValue":0.0},
                                               "m0":{"initialValue":0.0},
                                               "m1":{"initialValue":0.0}}

        Property = self.masstran.input.Property 
        Material = self.masstran.input.Material 

        i = 0
        for key, val in Property.items():
            i = i + 1
            Design_Variable_Relation = self.masstran.input.Design_Variable_Relation
            if Design_Variable_Relation is None: Design_Variable_Relation = {}
            
            if val["propertyType"] == "Shell":
                desvarR = {"componentType": "Property",
                           "componentName": key,
                           "fieldName" : "membraneThickness",
                           "constantCoeff" : val["membraneThickness"],
                           "variableName" : ["x0", "x1"],
                           "linearCoeff" : [2.0*i, 3.0*i]}

                Design_Variable_Relation["thick_"+key] = desvarR
            else:
                desvarR = {"componentType": "Property",
                           "componentName": key,
                           "fieldName" : "mass",
                           "constantCoeff" : val["mass"],
                           "variableName" : ["x0", "x1"],
                           "linearCoeff" : [1.1*i, 2.2*i]}

                Design_Variable_Relation["mass_"+key] = desvarR

            self.masstran.input.Design_Variable_Relation = Design_Variable_Relation

        i = 0
        for key, val in Material.items():
            i = i + 1
            Design_Variable_Relation = self.masstran.input.Design_Variable_Relation
            if Design_Variable_Relation is None: Design_Variable_Relation = {}
            
            desvarR = {"componentType": "Material",
                       "componentName": key,
                       "fieldName" : "density",
                       "constantCoeff" : val["density"],
                       "variableName" : ["m0", "m1"],
                       "linearCoeff" : [3.0*i, 4.0*i]}
    
            Design_Variable_Relation["rho_"+key] = desvarR
            self.masstran.input.Design_Variable_Relation = Design_Variable_Relation

        # derivatives of 1 box w.r.t pmtr
        Area_pmtr     = self.masstran.output["Area"].deriv()
        Mass_pmtr     = self.masstran.output["Mass"].deriv()
        Centroid_pmtr = self.masstran.output["Centroid"].deriv()
        CG_pmtr       = self.masstran.output["CG"].deriv()
        I_pmtr        = self.masstran.output["I_Vector"].deriv()
        IU_pmtr       = self.masstran.output["I_Upper"].deriv()
        IL_pmtr       = self.masstran.output["I_Lower"].deriv()
        II_pmtr       = self.masstran.output["I_Tensor"].deriv()

        self.masstran.input.Design_Variable = None

        for pmtr in ["L", "W", "H"]:
            print("===>", pmtr)
            # perturbed despmtr
            eps = 1e-4
            self.myProblem.geometry.despmtr[pmtr].value += eps
    
            AreaP     = self.masstran.output.Area
            MassP     = self.masstran.output.Mass
            CentroidP = self.masstran.output.Centroid
            CGP       = self.masstran.output.CG
            IP        = self.masstran.output.I_Vector
            IUP       = self.masstran.output.I_Upper
            ILP       = self.masstran.output.I_Lower
            IIP       = self.masstran.output.I_Tensor
            
            self.myProblem.geometry.despmtr[pmtr].value -= 2*eps
            
            AreaM     = self.masstran.output.Area
            MassM     = self.masstran.output.Mass
            CentroidM = self.masstran.output.Centroid
            CGM       = self.masstran.output.CG
            IM        = self.masstran.output.I_Vector
            IUM       = self.masstran.output.I_Upper
            ILM       = self.masstran.output.I_Lower
            IIM       = self.masstran.output.I_Tensor

            self.myProblem.geometry.despmtr[pmtr].value += eps

            self.assertAlmostEqual((AreaP-AreaM)/(2*eps), Area_pmtr[pmtr], 6)
            self.assertAlmostEqual((MassP-MassM)/(2*eps), Mass_pmtr[pmtr], 6)
            for i in range(len(CentroidP)):
                self.assertAlmostEqual((CentroidP[i]-CentroidM[i])/(2*eps), Centroid_pmtr[pmtr][i], 6)
            for i in range(len(CGP)):
                self.assertAlmostEqual((CGP[i]-CGM[i])/(2*eps), CG_pmtr[pmtr][i], 6)
            for i in range(len(IP)):
                self.assertAlmostEqual((IP[i]-IM[i])/(2*eps), I_pmtr[pmtr][i], Itol)
            for i in range(len(IUP)):
                self.assertAlmostEqual((IUP[i]-IUM[i])/(2*eps), IU_pmtr[pmtr][i], Itol)
            for i in range(len(ILP)):
                self.assertAlmostEqual((ILP[i]-ILM[i])/(2*eps), IL_pmtr[pmtr][i], Itol)
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual((IIP[i][j]-IIM[i][j])/(2*eps), II_pmtr[pmtr][3*i+j], Itol)


        for pmtr in ["x0", "x1", "m0", "m1"]:
            print("===>", pmtr)
            # perturbed Design_Varible
            eps = 1e-4
            
            Design_Variable = {"x0":{"initialValue":0.0}, "x1":{"initialValue":0.0},
                               "m0":{"initialValue":0.0}, "m1":{"initialValue":0.0}}
            
            Design_Variable[pmtr]["initialValue"] += eps
            self.masstran.input.Design_Variable = Design_Variable

            AreaP     = self.masstran.output.Area
            MassP     = self.masstran.output.Mass
            CentroidP = self.masstran.output.Centroid
            CGP       = self.masstran.output.CG
            IP        = self.masstran.output.I_Vector
            IUP       = self.masstran.output.I_Upper
            ILP       = self.masstran.output.I_Lower
            IIP       = self.masstran.output.I_Tensor
            
            Design_Variable[pmtr]["initialValue"] -= 2*eps
            self.masstran.input.Design_Variable = Design_Variable

            AreaM     = self.masstran.output.Area
            MassM     = self.masstran.output.Mass
            CentroidM = self.masstran.output.Centroid
            CGM       = self.masstran.output.CG
            IM        = self.masstran.output.I_Vector
            IUM       = self.masstran.output.I_Upper
            ILM       = self.masstran.output.I_Lower
            IIM       = self.masstran.output.I_Tensor

            self.assertAlmostEqual((AreaP-AreaM)/(2*eps), Area_pmtr[pmtr], 6)
            self.assertAlmostEqual((MassP-MassM)/(2*eps), Mass_pmtr[pmtr], 6)
            for i in range(len(CentroidP)):
                self.assertAlmostEqual((CentroidP[i]-CentroidM[i])/(2*eps), Centroid_pmtr[pmtr][i], 6)
            for i in range(len(CGP)):
                self.assertAlmostEqual((CGP[i]-CGM[i])/(2*eps), CG_pmtr[pmtr][i], 6)
            for i in range(len(IP)):
                self.assertAlmostEqual((IP[i]-IM[i])/(2*eps), I_pmtr[pmtr][i], 5)
            for i in range(len(IUP)):
                self.assertAlmostEqual((IUP[i]-IUM[i])/(2*eps), IU_pmtr[pmtr][i], 5)
            for i in range(len(ILP)):
                self.assertAlmostEqual((ILP[i]-ILM[i])/(2*eps), IL_pmtr[pmtr][i], 5)
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual((IIP[i][j]-IIM[i][j])/(2*eps), II_pmtr[pmtr][3*i+j], 5)

        self.masstran.input.Design_Variable = None
        self.masstran.input.Design_Variable_Relation = None

#==============================================================================
    def test_plate(self):

        # set the propery for this analysis
        self.masstran.input.Property = {"box1" : self.shell1}

        # Only work with plate
        self.myProblem.geometry.cfgpmtr["VIEW:plate"].value = 1
        self.myProblem.geometry.cfgpmtr["VIEW:box1" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box2" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box3" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:point"].value = 0

        # properties of 1 plate
        Area     = self.masstran.output.Area
        Mass     = self.masstran.output.Mass
        Centroid = self.masstran.output.Centroid
        CG       = self.masstran.output.CG
        I        = self.masstran.output.I_Vector
        IU       = self.masstran.output.I_Upper
        IL       = self.masstran.output.I_Lower
        II       = self.masstran.output.I_Tensor
        MassProp = self.masstran.output.MassProp

        self.assertEqual(3, len(Centroid))
        self.assertEqual(3, len(CG))
        self.assertEqual(6, len(I))
        self.assertEqual(3, len(II))
        for row in II:
            self.assertEqual(3, len(row))

        Ixx = I[0]; Iyy = I[1]; Izz = I[2]
        Ixy = I[3]; Ixz = I[4]; Iyz = I[5]

        areaTrue = 2*3
        m = 2*10*areaTrue # Mass of the plate
        IxxTrue = 1./12.*m*3**2
        IyyTrue = 1./12.*m*2**2
        IzzTrue = 1./12.*m*(2**2 + 3**2)

        self.assertAlmostEqual(areaTrue, Area, 9)
        self.assertAlmostEqual(m, Mass, 9)

        self.assertAlmostEqual(1.0, Centroid[0], 9)
        self.assertAlmostEqual(1.5, Centroid[1], 9)
        self.assertAlmostEqual(0.0, Centroid[2], 9)

        self.assertAlmostEqual(1.0, CG[0], 9)
        self.assertAlmostEqual(1.5, CG[1], 9)
        self.assertAlmostEqual(0.0, CG[2], 9)

        self.assertAlmostEqual(Ixx/IxxTrue, 1.0, 2)
        self.assertAlmostEqual(Iyy/IyyTrue, 1.0, 2)
        self.assertAlmostEqual(Izz/IzzTrue, 1.0, 2)

        self.assertAlmostEqual(0.0, Ixy, 9)
        self.assertAlmostEqual(0.0, Ixz, 9)
        self.assertAlmostEqual(0.0, Iyz, 9)

        self.assertAlmostEqual(Ixx, self.masstran.output.Ixx, 9)
        self.assertAlmostEqual(Iyy, self.masstran.output.Iyy, 9)
        self.assertAlmostEqual(Izz, self.masstran.output.Izz, 9)
        self.assertAlmostEqual(Ixy, self.masstran.output.Ixy, 9)
        self.assertAlmostEqual(Ixz, self.masstran.output.Ixz, 9)
        self.assertAlmostEqual(Iyz, self.masstran.output.Iyz, 9)

        self.assertAlmostEqual( Ixx, IU[0], 9)
        self.assertAlmostEqual(-Ixy, IU[1], 9)
        self.assertAlmostEqual(-Ixz, IU[2], 9)
        self.assertAlmostEqual( Iyy, IU[3], 9)
        self.assertAlmostEqual(-Iyz, IU[4], 9)
        self.assertAlmostEqual( Izz, IU[5], 9)

        self.assertAlmostEqual( Ixx, IL[0], 9)
        self.assertAlmostEqual(-Ixy, IL[1], 9)
        self.assertAlmostEqual( Iyy, IL[2], 9)
        self.assertAlmostEqual(-Ixz, IL[3], 9)
        self.assertAlmostEqual(-Iyz, IL[4], 9)
        self.assertAlmostEqual( Izz, IL[5], 9)

        self.assertAlmostEqual( Ixx, II[0][0], 9)
        self.assertAlmostEqual(-Ixy, II[0][1], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Ixy, II[1][0], 9)
        self.assertAlmostEqual( Iyy, II[1][1], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual( Izz, II[2][2], 9)

        #self.assertAlmostEqual(m  , MassProp["mass"],  9)
        #self.assertAlmostEqual(1.0, MassProp["CG"][0], 9)
        #self.assertAlmostEqual(1.5, MassProp["CG"][1], 9)
        #self.assertAlmostEqual(0.0, MassProp["CG"][2], 9)
        #self.assertAlmostEqual(Ixx, MassProp["massInertia"][0], 9)
        #self.assertAlmostEqual(Iyy, MassProp["massInertia"][1], 9)
        #self.assertAlmostEqual(Izz, MassProp["massInertia"][2], 9)
        #self.assertAlmostEqual(Ixy, MassProp["massInertia"][3], 9)
        #self.assertAlmostEqual(Ixz, MassProp["massInertia"][4], 9)
        #self.assertAlmostEqual(Iyz, MassProp["massInertia"][5], 9)
        
        self.check_derivs()


#==============================================================================
    def test_plate_point(self):

        # set the propery for this analysis
        self.masstran.input.Property = {"box1" : self.shell1, "point" : self.pointMass}

        despmtr = self.myProblem.geometry.despmtr

        # Only work with plate
        self.myProblem.geometry.cfgpmtr["VIEW:plate"].value = 1
        self.myProblem.geometry.cfgpmtr["VIEW:box1" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box2" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box3" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:point"].value = 1

        # properties of 1 plate
        Area     = self.masstran.output.Area
        Mass     = self.masstran.output.Mass
        Centroid = self.masstran.output.Centroid
        CG       = self.masstran.output.CG
        I        = self.masstran.output.I_Vector
        IU       = self.masstran.output.I_Upper
        IL       = self.masstran.output.I_Lower
        II       = self.masstran.output.I_Tensor
        MassProp = self.masstran.output.MassProp

        self.assertEqual(3, len(Centroid))
        self.assertEqual(3, len(CG))
        self.assertEqual(6, len(I))
        self.assertEqual(3, len(II))
        for row in II:
            self.assertEqual(3, len(row))

        Ixx = I[0]; Iyy = I[1]; Izz = I[2]
        Ixy = I[3]; Ixz = I[4]; Iyz = I[5]

        areaTrue = 2.0*despmtr.L*3.0*despmtr.W
        mPlate  = 2*10*areaTrue # Mass of the plate
        
        xcentPl = 2.0*despmtr.L/2
        ycentPl = 3.0*despmtr.W/2
        zcentPl = 0.0

        # Add plate inertia        
        IxxTrue = 4./12.*mPlate*3**2
        IyyTrue = 4./12.*mPlate*2**2
        IzzTrue = 4./12.*mPlate*(2**2 + 3**2)
        IxyTrue = 3./12.*mPlate*(2*3)
        IxzTrue = 0
        IyzTrue = 0
       
        mPoint = self.pointMass["mass"] # add point mass
        
        m = mPlate + mPoint
        
        # Point centroid
        xcentPn = despmtr.L
        ycentPn = despmtr.W
        zcentPn = despmtr.H
        
        # Overall CG
        CGTrue = [(xcentPl*mPlate + xcentPn*mPoint)/m, 
                  (ycentPl*mPlate + ycentPn*mPoint)/m, 
                  (zcentPl*mPlate + zcentPn*mPoint)/m]

        massInertia = self.pointMass["massInertia"]

        # Add point inertia
        IxxTrue += (ycentPn**2 + zcentPn**2) * mPoint + massInertia[0]
        IyyTrue += (xcentPn**2 + zcentPn**2) * mPoint + massInertia[2]
        IzzTrue += (xcentPn**2 + ycentPn**2) * mPoint + massInertia[5]
        IxyTrue += (xcentPn    * ycentPn   ) * mPoint - massInertia[1]
        IxzTrue += (xcentPn    * zcentPn   ) * mPoint - massInertia[3]
        IyzTrue += (ycentPn    * zcentPn   ) * mPoint - massInertia[4]

        # overall inertia
        IxxTrue += - m * (CGTrue[1]**2 + CGTrue[2]**2)
        IyyTrue += - m * (CGTrue[0]**2 + CGTrue[2]**2)
        IzzTrue += - m * (CGTrue[0]**2 + CGTrue[1]**2)
        IxyTrue += - m *  CGTrue[0] * CGTrue[1]
        IxzTrue += - m *  CGTrue[0] * CGTrue[2]
        IyzTrue += - m *  CGTrue[1] * CGTrue[2]

        self.assertAlmostEqual(areaTrue, Area, 9)
        self.assertAlmostEqual(m, Mass, 9)

        self.assertAlmostEqual(1.0, Centroid[0], 9)
        self.assertAlmostEqual(1.5, Centroid[1], 9)
        self.assertAlmostEqual(0.0, Centroid[2], 9)

        self.assertAlmostEqual(CGTrue[0], CG[0], 9)
        self.assertAlmostEqual(CGTrue[1], CG[1], 9)
        self.assertAlmostEqual(CGTrue[2], CG[2], 9)

        self.assertAlmostEqual(Ixx/IxxTrue, 1.0, 2)
        self.assertAlmostEqual(Iyy/IyyTrue, 1.0, 2)
        self.assertAlmostEqual(Izz/IzzTrue, 1.0, 2)

        self.assertAlmostEqual(Ixy/IxyTrue, 1.0, 9)
        self.assertAlmostEqual(Ixz/IxzTrue, 1.0, 9)
        self.assertAlmostEqual(Iyz/IyzTrue, 1.0, 9)

        self.assertAlmostEqual(Ixx, self.masstran.output.Ixx, 9)
        self.assertAlmostEqual(Iyy, self.masstran.output.Iyy, 9)
        self.assertAlmostEqual(Izz, self.masstran.output.Izz, 9)
        self.assertAlmostEqual(Ixy, self.masstran.output.Ixy, 9)
        self.assertAlmostEqual(Ixz, self.masstran.output.Ixz, 9)
        self.assertAlmostEqual(Iyz, self.masstran.output.Iyz, 9)

        self.assertAlmostEqual( Ixx, IU[0], 9)
        self.assertAlmostEqual(-Ixy, IU[1], 9)
        self.assertAlmostEqual(-Ixz, IU[2], 9)
        self.assertAlmostEqual( Iyy, IU[3], 9)
        self.assertAlmostEqual(-Iyz, IU[4], 9)
        self.assertAlmostEqual( Izz, IU[5], 9)

        self.assertAlmostEqual( Ixx, IL[0], 9)
        self.assertAlmostEqual(-Ixy, IL[1], 9)
        self.assertAlmostEqual( Iyy, IL[2], 9)
        self.assertAlmostEqual(-Ixz, IL[3], 9)
        self.assertAlmostEqual(-Iyz, IL[4], 9)
        self.assertAlmostEqual( Izz, IL[5], 9)

        self.assertAlmostEqual( Ixx, II[0][0], 9)
        self.assertAlmostEqual(-Ixy, II[0][1], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Ixy, II[1][0], 9)
        self.assertAlmostEqual( Iyy, II[1][1], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual( Izz, II[2][2], 9)

        #self.assertAlmostEqual(m  , MassProp["mass"],  9)
        #self.assertAlmostEqual(1.0, MassProp["CG"][0], 9)
        #self.assertAlmostEqual(1.5, MassProp["CG"][1], 9)
        #self.assertAlmostEqual(0.0, MassProp["CG"][2], 9)
        #self.assertAlmostEqual(Ixx, MassProp["massInertia"][0], 9)
        #self.assertAlmostEqual(Iyy, MassProp["massInertia"][1], 9)
        #self.assertAlmostEqual(Izz, MassProp["massInertia"][2], 9)
        #self.assertAlmostEqual(Ixy, MassProp["massInertia"][3], 9)
        #self.assertAlmostEqual(Ixz, MassProp["massInertia"][4], 9)
        #self.assertAlmostEqual(Iyz, MassProp["massInertia"][5], 9)
        
        self.check_derivs()

#==============================================================================
    def test_box1(self):

        # set the propery for this analysis
        self.masstran.input.Property = {"box1" : self.shell1}

        # Only work with 1 box now
        self.myProblem.geometry.cfgpmtr["VIEW:plate"].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box1" ].value = 1
        self.myProblem.geometry.cfgpmtr["VIEW:box2" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box3" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:point"].value = 0

        # properties of 1 box
        Area     = self.masstran.output.Area
        Mass     = self.masstran.output.Mass
        Centroid = self.masstran.output.Centroid
        CG       = self.masstran.output.CG
        I        = self.masstran.output.I_Vector
        IU       = self.masstran.output.I_Upper
        IL       = self.masstran.output.I_Lower
        II       = self.masstran.output.I_Tensor
        MassProp = self.masstran.output.MassProp

        self.assertEqual(3, len(Centroid))
        self.assertEqual(3, len(CG))
        self.assertEqual(6, len(I))
        self.assertEqual(3, len(II))
        for row in II:
            self.assertEqual(3, len(row))

        Ixx = I[0]; Iyy = I[1]; Izz = I[2]
        Ixy = I[3]; Ixz = I[4]; Iyz = I[5]

        nFace = 6
        areaTrue = 1
        m = 2*10 # Mass of one face of the box
        massTrue = m*nFace
        
        # CG location of each face
        fCG = {0:[0.0, 0.5, 0.5], 1:[1.0, 0.5, 0.5], # x-min/max
               2:[0.5, 0.0, 0.5], 3:[0.5, 1.0, 0.5], # y-min/max
               4:[0.5, 0.5, 0.0], 5:[0.5, 0.5, 1.0]} # z-min/max

        # length of each face
        fLen = {0:[0.0, 1.0, 1.0], 1:[0.0, 1.0, 1.0], # x-min/max
                2:[1.0, 0.0, 1.0], 3:[1.0, 0.0, 1.0], # y-min/max
                4:[1.0, 1.0, 0.0], 5:[1.0, 1.0, 0.0]} # z-min/max

        IxxTrue = 0
        IyyTrue = 0
        IzzTrue = 0
        for i in range(6):
            IxxTrue += m*(fCG[i][1]**2 + fCG[i][2]**2) + 1/12.*m*(fLen[i][1]**2 + fLen[i][2]**2)
            IyyTrue += m*(fCG[i][0]**2 + fCG[i][2]**2) + 1/12.*m*(fLen[i][0]**2 + fLen[i][2]**2)
            IzzTrue += m*(fCG[i][0]**2 + fCG[i][1]**2) + 1/12.*m*(fLen[i][0]**2 + fLen[i][1]**2)

        IxxTrue -= massTrue*(CG[1]**2 + CG[2]**2)
        IyyTrue -= massTrue*(CG[0]**2 + CG[2]**2)
        IzzTrue -= massTrue*(CG[0]**2 + CG[1]**2)

        self.assertAlmostEqual(nFace*areaTrue, Area, 9)
        self.assertAlmostEqual(nFace*m, Mass, 9)

        self.assertAlmostEqual(0.5, Centroid[0], 9)
        self.assertAlmostEqual(0.5, Centroid[1], 9)
        self.assertAlmostEqual(0.5, Centroid[2], 9)

        self.assertAlmostEqual(0.5, CG[0], 9)
        self.assertAlmostEqual(0.5, CG[1], 9)
        self.assertAlmostEqual(0.5, CG[2], 9)

        self.assertAlmostEqual(Ixx/IxxTrue, 1.0, 2)
        self.assertAlmostEqual(Iyy/IyyTrue, 1.0, 2)
        self.assertAlmostEqual(Izz/IzzTrue, 1.0, 2)

        self.assertAlmostEqual(0.0, Ixy, 9)
        self.assertAlmostEqual(0.0, Ixz, 9)
        self.assertAlmostEqual(0.0, Iyz, 9)

        self.assertAlmostEqual(Ixx, self.masstran.output.Ixx, 9)
        self.assertAlmostEqual(Iyy, self.masstran.output.Iyy, 9)
        self.assertAlmostEqual(Izz, self.masstran.output.Izz, 9)
        self.assertAlmostEqual(Ixy, self.masstran.output.Ixy, 9)
        self.assertAlmostEqual(Ixz, self.masstran.output.Ixz, 9)
        self.assertAlmostEqual(Iyz, self.masstran.output.Iyz, 9)

        self.assertAlmostEqual( Ixx, IU[0], 9)
        self.assertAlmostEqual(-Ixy, IU[1], 9)
        self.assertAlmostEqual(-Ixz, IU[2], 9)
        self.assertAlmostEqual( Iyy, IU[3], 9)
        self.assertAlmostEqual(-Iyz, IU[4], 9)
        self.assertAlmostEqual( Izz, IU[5], 9)

        self.assertAlmostEqual( Ixx, IL[0], 9)
        self.assertAlmostEqual(-Ixy, IL[1], 9)
        self.assertAlmostEqual( Iyy, IL[2], 9)
        self.assertAlmostEqual(-Ixz, IL[3], 9)
        self.assertAlmostEqual(-Iyz, IL[4], 9)
        self.assertAlmostEqual( Izz, IL[5], 9)

        self.assertAlmostEqual( Ixx, II[0][0], 9)
        self.assertAlmostEqual(-Ixy, II[0][1], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Ixy, II[1][0], 9)
        self.assertAlmostEqual( Iyy, II[1][1], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual( Izz, II[2][2], 9)

        #self.assertAlmostEqual(m  , MassProp["mass"],  9)
        #self.assertAlmostEqual(0.5, MassProp["CG"][0], 9)
        #self.assertAlmostEqual(0.5, MassProp["CG"][1], 9)
        #self.assertAlmostEqual(0.5, MassProp["CG"][2], 9)
        #self.assertAlmostEqual(Ixx, MassProp["massInertia"][0], 9)
        #self.assertAlmostEqual(Iyy, MassProp["massInertia"][1], 9)
        #self.assertAlmostEqual(Izz, MassProp["massInertia"][2], 9)
        #self.assertAlmostEqual(Ixy, MassProp["massInertia"][3], 9)
        #self.assertAlmostEqual(Ixz, MassProp["massInertia"][4], 9)
        #self.assertAlmostEqual(Iyz, MassProp["massInertia"][5], 9)

        self.check_derivs()

#==============================================================================
    def test_box2(self):
 
 
        egadsTess = self.myProblem.analysis.create(aim = "egadsTessAIM",
                                                   name = "egadsTess")

        egadsTess.input.Mesh_Quiet_Flag = True
        egadsTess.input.Tess_Params = [0.05, 0.01, 15.0]

        self.masstran.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])
        self.masstran.input.Mesh_Morph = True

        # set the propery for this analysis
        self.masstran.input.Property = {"box1" : self.shell1}

        # Only work with 2 box now
        self.myProblem.geometry.cfgpmtr["VIEW:plate"].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box1" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box2" ].value = 1
        self.myProblem.geometry.cfgpmtr["VIEW:box3" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:point"].value = 0

        # properties of 1 box lacking 1 face
        Area     = self.masstran.output.Area
        Mass     = self.masstran.output.Mass
        Centroid = self.masstran.output.Centroid
        CG       = self.masstran.output.CG
        I        = self.masstran.output.I_Vector
        IU       = self.masstran.output.I_Upper
        IL       = self.masstran.output.I_Lower
        II       = self.masstran.output.I_Tensor
        MassProp = self.masstran.output.MassProp
 
        self.assertEqual(3, len(Centroid))
        self.assertEqual(3, len(CG))
        self.assertEqual(6, len(I))
        self.assertEqual(3, len(II))
        for row in II:
            self.assertEqual(3, len(row))

        Ixx = I[0]; Iyy = I[1]; Izz = I[2]
        Ixy = I[3]; Ixz = I[4]; Iyz = I[5]

        nFace = 6
        areaTrue = 1
        m = 2*10 # Mass of one face of the box
        massTrue = m*nFace
        
        # CG location of each face
        fCG = {0:[2.0, 0.5, 0.5], 1:[3.0, 0.5, 0.5], # x-min/max
               2:[2.5, 0.0, 0.5], 3:[2.5, 1.0, 0.5], # y-min/max
               4:[2.5, 0.5, 0.0], 5:[2.5, 0.5, 1.0]} # z-min/max

        # length of each face
        fLen = {0:[0.0, 1.0, 1.0], 1:[0.0, 1.0, 1.0], # x-min/max
                2:[1.0, 0.0, 1.0], 3:[1.0, 0.0, 1.0], # y-min/max
                4:[1.0, 1.0, 0.0], 5:[1.0, 1.0, 0.0]} # z-min/max

        IxxTrue = 0
        IyyTrue = 0
        IzzTrue = 0
        IxyTrue = 0
        IxzTrue = 0
        IyzTrue = 0
        for i in range(6-nFace,6):
            IxxTrue += m*(fCG[i][1]**2 + fCG[i][2]**2) + 1/12.*m*(fLen[i][1]**2 + fLen[i][2]**2)
            IyyTrue += m*(fCG[i][0]**2 + fCG[i][2]**2) + 1/12.*m*(fLen[i][0]**2 + fLen[i][2]**2)
            IzzTrue += m*(fCG[i][0]**2 + fCG[i][1]**2) + 1/12.*m*(fLen[i][0]**2 + fLen[i][1]**2)

            IxyTrue += m*(fCG[i][0]*fCG[i][1])
            IxzTrue += m*(fCG[i][0]*fCG[i][2])
            IyzTrue += m*(fCG[i][1]*fCG[i][2])

        IxxTrue -= massTrue*(CG[1]**2 + CG[2]**2)
        IyyTrue -= massTrue*(CG[0]**2 + CG[2]**2)
        IzzTrue -= massTrue*(CG[0]**2 + CG[1]**2)

        IxyTrue -= massTrue*(CG[0]*CG[1])
        IxzTrue -= massTrue*(CG[0]*CG[2])
        IyzTrue -= massTrue*(CG[1]*CG[2])
        
        self.assertAlmostEqual(nFace*areaTrue, Area, 9)
        self.assertAlmostEqual(massTrue, Mass, 9)
 
        # compute centroid as averge of last 5 faces plus the shift by 2
        Cx = 2+(4*0.5+1)/nFace
 
        self.assertAlmostEqual(Cx , Centroid[0], 9)
        self.assertAlmostEqual(0.5, Centroid[1], 9)
        self.assertAlmostEqual(0.5, Centroid[2], 9)
 
        self.assertAlmostEqual(Cx,  CG[0], 9)
        self.assertAlmostEqual(0.5, CG[1], 9)
        self.assertAlmostEqual(0.5, CG[2], 9)

        self.assertAlmostEqual(Ixx/IxxTrue, 1.0, 2)
        self.assertAlmostEqual(Iyy/IyyTrue, 1.0, 2)
        self.assertAlmostEqual(Izz/IzzTrue, 1.0, 2)

        self.assertAlmostEqual(IxyTrue, Ixy, 2)
        self.assertAlmostEqual(IxzTrue, Ixz, 2)
        self.assertAlmostEqual(IyzTrue, Iyz, 2)

        self.assertAlmostEqual(Ixx, self.masstran.output.Ixx, 9)
        self.assertAlmostEqual(Iyy, self.masstran.output.Iyy, 9)
        self.assertAlmostEqual(Izz, self.masstran.output.Izz, 9)
        self.assertAlmostEqual(Ixy, self.masstran.output.Ixy, 9)
        self.assertAlmostEqual(Ixz, self.masstran.output.Ixz, 9)
        self.assertAlmostEqual(Iyz, self.masstran.output.Iyz, 9)

        self.assertAlmostEqual( Ixx, IU[0], 9)
        self.assertAlmostEqual(-Ixy, IU[1], 9)
        self.assertAlmostEqual(-Ixz, IU[2], 9)
        self.assertAlmostEqual( Iyy, IU[3], 9)
        self.assertAlmostEqual(-Iyz, IU[4], 9)
        self.assertAlmostEqual( Izz, IU[5], 9)

        self.assertAlmostEqual( Ixx, IL[0], 9)
        self.assertAlmostEqual(-Ixy, IL[1], 9)
        self.assertAlmostEqual( Iyy, IL[2], 9)
        self.assertAlmostEqual(-Ixz, IL[3], 9)
        self.assertAlmostEqual(-Iyz, IL[4], 9)
        self.assertAlmostEqual( Izz, IL[5], 9)

        self.assertAlmostEqual( Ixx, II[0][0], 9)
        self.assertAlmostEqual(-Ixy, II[0][1], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Ixy, II[1][0], 9)
        self.assertAlmostEqual( Iyy, II[1][1], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual( Izz, II[2][2], 9)

        #self.assertAlmostEqual(m  , MassProp["mass"],  9)
        #self.assertAlmostEqual(Cx , MassProp["CG"][0], 9)
        #self.assertAlmostEqual(0.5, MassProp["CG"][1], 9)
        #self.assertAlmostEqual(0.5, MassProp["CG"][2], 9)
        #self.assertAlmostEqual(Ixx, MassProp["massInertia"][0], 9)
        #self.assertAlmostEqual(Iyy, MassProp["massInertia"][1], 9)
        #self.assertAlmostEqual(Izz, MassProp["massInertia"][2], 9)
        #self.assertAlmostEqual(Ixy, MassProp["massInertia"][3], 9)
        #self.assertAlmostEqual(Ixz, MassProp["massInertia"][4], 9)
        #self.assertAlmostEqual(Iyz, MassProp["massInertia"][5], 9)

        self.masstran.input["Surface_Mesh"].unlink()
        self.check_derivs(Itol=2)
        self.masstran.input.Mesh_Morph = False

#==============================================================================
    def test_box3(self):
 
        # set the propery for this analysis
        self.masstran.input.Property = {"box1" : self.shell1, "box2" : self.shell2}

        # Only work with 3 box now
        self.myProblem.geometry.cfgpmtr["VIEW:plate"].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box1" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box2" ].value = 0
        self.myProblem.geometry.cfgpmtr["VIEW:box3" ].value = 1
        self.myProblem.geometry.cfgpmtr["VIEW:point"].value = 0

        # properties of 1 box lacking 1 face
        Area     = self.masstran.output.Area
        Mass     = self.masstran.output.Mass
        Centroid = self.masstran.output.Centroid
        CG       = self.masstran.output.CG
        I        = self.masstran.output.I_Vector
        II       = self.masstran.output.I_Tensor
        MassProp = self.masstran.output.MassProp

        self.assertEqual(3, len(Centroid))
        self.assertEqual(3, len(CG))
        self.assertEqual(6, len(I))
        for row in II:
            self.assertEqual(3, len(row))

        Ixx = I[0]; Iyy = I[1]; Izz = I[2]
        Ixy = I[3]; Ixz = I[4]; Iyz = I[5]

        nFace = 6
        areaTrue = 1
        
        # Mass of faces of the box
        m = {0:2*10, 1:3*20,
             2:2*10, 3:3*20,
             4:2*10, 5:3*20}
        
        # total Mass
        massTrue = 0
        for i in range(6):
            massTrue += m[i]
        
        # CG location of each face
        fCG = {0:[4.0, 0.5, 0.5], 1:[5.0, 0.5, 0.5], # x-min/max
               2:[4.5, 0.0, 0.5], 3:[4.5, 1.0, 0.5], # y-min/max
               4:[4.5, 0.5, 0.0], 5:[4.5, 0.5, 1.0]} # z-min/max

        # length of each face
        fLen = {0:[0.0, 1.0, 1.0], 1:[0.0, 1.0, 1.0], # x-min/max
                2:[1.0, 0.0, 1.0], 3:[1.0, 0.0, 1.0], # y-min/max
                4:[1.0, 1.0, 0.0], 5:[1.0, 1.0, 0.0]} # z-min/max

        IxxTrue = 0
        IyyTrue = 0
        IzzTrue = 0
        IxyTrue = 0
        IxzTrue = 0
        IyzTrue = 0
        for i in range(6):
            IxxTrue += m[i]*(fCG[i][1]**2 + fCG[i][2]**2) + 1/12.*m[i]*(fLen[i][1]**2 + fLen[i][2]**2)
            IyyTrue += m[i]*(fCG[i][0]**2 + fCG[i][2]**2) + 1/12.*m[i]*(fLen[i][0]**2 + fLen[i][2]**2)
            IzzTrue += m[i]*(fCG[i][0]**2 + fCG[i][1]**2) + 1/12.*m[i]*(fLen[i][0]**2 + fLen[i][1]**2)

            IxyTrue += m[i]*(fCG[i][0]*fCG[i][1])
            IxzTrue += m[i]*(fCG[i][0]*fCG[i][2])
            IyzTrue += m[i]*(fCG[i][1]*fCG[i][2])

        IxxTrue -= massTrue*(CG[1]**2 + CG[2]**2)
        IyyTrue -= massTrue*(CG[0]**2 + CG[2]**2)
        IzzTrue -= massTrue*(CG[0]**2 + CG[1]**2)

        IxyTrue -= massTrue*(CG[0]*CG[1])
        IxzTrue -= massTrue*(CG[0]*CG[2])
        IyzTrue -= massTrue*(CG[1]*CG[2])

        self.assertAlmostEqual(nFace*areaTrue, Area, 9)
        self.assertAlmostEqual(massTrue, Mass, 9)
 
        self.assertAlmostEqual(4.5, Centroid[0], 9)
        self.assertAlmostEqual(0.5, Centroid[1], 9)
        self.assertAlmostEqual(0.5, Centroid[2], 9)
  
        CGx = 0
        CGy = 0
        CGz = 0
        for i in range(6):
            CGx += fCG[i][0]*m[i]
            CGy += fCG[i][1]*m[i]
            CGz += fCG[i][2]*m[i]
        CGx /= massTrue
        CGy /= massTrue
        CGz /= massTrue

        self.assertAlmostEqual(CGx, CG[0], 9)
        self.assertAlmostEqual(CGy, CG[1], 9)
        self.assertAlmostEqual(CGz, CG[2], 9)

        self.assertAlmostEqual(Ixx/IxxTrue, 1.0, 2)
        self.assertAlmostEqual(Iyy/IyyTrue, 1.0, 2)
        self.assertAlmostEqual(Izz/IzzTrue, 1.0, 2)

        self.assertAlmostEqual(IxyTrue, Ixy, 4)
        self.assertAlmostEqual(IxzTrue, Ixz, 4)
        self.assertAlmostEqual(IyzTrue, Iyz, 4)

        self.assertAlmostEqual(Ixx, self.masstran.output.Ixx, 9)
        self.assertAlmostEqual(Iyy, self.masstran.output.Iyy, 9)
        self.assertAlmostEqual(Izz, self.masstran.output.Izz, 9)
        self.assertAlmostEqual(Ixy, self.masstran.output.Ixy, 9)
        self.assertAlmostEqual(Ixz, self.masstran.output.Ixz, 9)
        self.assertAlmostEqual(Iyz, self.masstran.output.Iyz, 9)

        self.assertAlmostEqual( Ixx, II[0][0], 9)
        self.assertAlmostEqual(-Ixy, II[0][1], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Ixy, II[1][0], 9)
        self.assertAlmostEqual( Iyy, II[1][1], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual(-Ixz, II[0][2], 9)
        self.assertAlmostEqual(-Iyz, II[1][2], 9)
        self.assertAlmostEqual( Izz, II[2][2], 9)

        #self.assertAlmostEqual(m  , MassProp["mass"],  9)
        #self.assertAlmostEqual(CGx, MassProp["CG"][0], 9)
        #self.assertAlmostEqual(CGy, MassProp["CG"][1], 9)
        #self.assertAlmostEqual(CGz, MassProp["CG"][2], 9)
        #self.assertAlmostEqual(Ixx, MassProp["massInertia"][0], 9)
        #self.assertAlmostEqual(Iyy, MassProp["massInertia"][1], 9)
        #self.assertAlmostEqual(Izz, MassProp["massInertia"][2], 9)
        #self.assertAlmostEqual(Ixy, MassProp["massInertia"][3], 9)
        #self.assertAlmostEqual(Ixz, MassProp["massInertia"][4], 9)
        #self.assertAlmostEqual(Iyz, MassProp["massInertia"][5], 9)

        self.check_derivs()

#==============================================================================
    def run_journal(self, myProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Only work with 2 boxes
        if verbose: print(6*"-","VIEW:plate", line)
        myProblem.geometry.cfgpmtr["VIEW:plate"].value = 0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","VIEW:box3", line)
        myProblem.geometry.cfgpmtr["VIEW:box3" ].value = 0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","VIEW:point", line)
        myProblem.geometry.cfgpmtr["VIEW:point"].value = 1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load masstran AIM
        if verbose: print(6*"-","Load masstranAIM", line)
        masstran = myProblem.analysis.create(aim = "masstranAIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Set edge points
        if verbose: print(6*"-","Modify Edge_Point_Min", line)
        masstran.input.Edge_Point_Min = 20; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Modify Edge_Point_Max", line)
        masstran.input.Edge_Point_Max = 20; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # set the Property for this analysis
        if verbose: print(6*"-","Modify Property", line)
        masstran.input.Property = {"box1" : self.shell1, "point" : self.pointMass}; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # set the Material for this analysis
        if verbose: print(6*"-","Modify Material", line)
        masstran.input.Material = self.masstran.input.Material; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # properties
        if verbose: print(6*"-","Area_1", line)
        Area_1 = masstran.output.Area; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Mass_1", line)
        Mass_1 = masstran.output.Mass; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","L", line)
        myProblem.geometry.despmtr.L *= 2; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

       # properties
        if verbose: print(6*"-","Area_2", line)
        Area_2 = masstran.output.Area; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Mass_2", line)
        Mass_2 = masstran.output.Mass; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())


        # Check that the counts have decreased
        self.assertGreater(Area_2, Area_1)
        self.assertGreater(Mass_2, Mass_1)

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def test_journal(self):

        capsFile = os.path.join("..","csmData","masstran.csm")
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
