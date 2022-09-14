# Import other need modules
from __future__ import print_function

import unittest
import os
import glob
import shutil

# Import pyCAPS class file
import pyCAPS


# Geometry is verified with avl by plotting the camber using the commands:
#
# avl caps
# oper
# g
# ca

class TestAVL(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # Create working directory variable
        cls.problemName = "workDir_avlAnalysisTest"

        cls.cleanUp()

        # Initialize Problem object
        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile="../csmData/avlSections.csm", outLevel=0)

    @classmethod
    def tearDownClass(cls):
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove problemName directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_numSpan(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanTotal"      : 24, # Can specify both Total and PerSection
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        wing = {"groupName"  : "Wing",
                "numChord"   : 8,
                "spaceChord" : 1.0,
                "numSpan"    : 24, # numSpan is depricated
                "spaceSpan"  : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    def test_alpha_custom_increment(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        Alpha      = [0.0, 3.0, 9.0]
        CLtotTrue  = [ 0.23858,  0.41257,  0.74994]
        CDtotTrue  = [ 0.01285,  0.02492,  0.06719]
        CmtotTrue  = [-0.10501, -0.09945, -0.08685]

        for i in range(0,len(Alpha)):
            # Set custom AoA
            avl.input.Alpha = Alpha[i]

            # Retrieve results
            CLtot = avl.output.CLtot
            CDtot = avl.output.CDtot
            Cmtot = avl.output.Cmtot
            #print("Alpha = ", Alpha[i])
            #print("CLtot = ", CLtot)
            #print("CDtot = ", CDtot)
            #print("Cmtot = ", Cmtot)

            self.assertAlmostEqual(CLtotTrue[i], CLtot, 4)
            self.assertAlmostEqual(CDtotTrue[i], CDtot, 4)
            self.assertAlmostEqual(CmtotTrue[i], Cmtot, 4)

            # Check derivatives

            # Stability axis angles
            self.assertAlmostEqual(avl.output["CLa"].value, avl.output["CLtot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["CLb"].value, avl.output["CLtot"].deriv("Beta"), 8)

            self.assertAlmostEqual(avl.output["Cl'a"].value, avl.output["Cl'tot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["Cl'b"].value, avl.output["Cl'tot"].deriv("Beta"), 8)

            self.assertAlmostEqual(avl.output["Cma"].value, avl.output["Cmtot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["Cmb"].value, avl.output["Cmtot"].deriv("Beta"), 8)

            self.assertAlmostEqual(avl.output["Cn'a"].value, avl.output["Cn'tot"].deriv("Alpha"), 8)
            self.assertAlmostEqual(avl.output["Cn'b"].value, avl.output["Cn'tot"].deriv("Beta"), 8)

            # Body axis rotation rates
            self.assertAlmostEqual(avl.output["CXp"].value, avl.output["CXtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["CXq"].value, avl.output["CXtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["CXr"].value, avl.output["CXtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["CYp"].value, avl.output["CYtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["CYq"].value, avl.output["CYtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["CYr"].value, avl.output["CYtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["CZp"].value, avl.output["CZtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["CZq"].value, avl.output["CZtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["CZr"].value, avl.output["CZtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["Clp"].value, avl.output["Cltot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["Clq"].value, avl.output["Cltot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["Clr"].value, avl.output["Cltot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["Cmp"].value, avl.output["Cmtot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["Cmq"].value, avl.output["Cmtot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["Cmr"].value, avl.output["Cmtot"].deriv("YawRate"), 8)

            self.assertAlmostEqual(avl.output["Cnp"].value, avl.output["Cntot"].deriv("RollRate"), 8)
            self.assertAlmostEqual(avl.output["Cnq"].value, avl.output["Cntot"].deriv("PitchRate"), 8)
            self.assertAlmostEqual(avl.output["Cnr"].value, avl.output["Cntot"].deriv("YawRate"), 8)

            # Control surfaces in the stability axis
            ControlStability = avl.output["ControlStability"].value

            self.assertAlmostEqual(ControlStability['LeftAileron']["CLtot"] , avl.output["CLtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["CYtot"] , avl.output["CYtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["Cl'tot"], avl.output["Cl'tot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["Cmtot"] , avl.output["Cmtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlStability['LeftAileron']["Cn'tot"], avl.output["Cn'tot"].deriv("LeftAileron"), 8)

            self.assertAlmostEqual(ControlStability['RightAileron']["CLtot"] , avl.output["CLtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["CYtot"] , avl.output["CYtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["Cl'tot"], avl.output["Cl'tot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["Cmtot"] , avl.output["Cmtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlStability['RightAileron']["Cn'tot"], avl.output["Cn'tot"].deriv("RightAileron"), 8)


            # Control surfaces in the body axis
            ControlBody = avl.output["ControlBody"].value

            self.assertAlmostEqual(ControlBody['LeftAileron']["CXtot"], avl.output["CXtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["CYtot"], avl.output["CYtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["CZtot"], avl.output["CZtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["Cltot"], avl.output["Cltot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["Cmtot"], avl.output["Cmtot"].deriv("LeftAileron"), 8)
            self.assertAlmostEqual(ControlBody['LeftAileron']["Cntot"], avl.output["Cntot"].deriv("LeftAileron"), 8)

            self.assertAlmostEqual(ControlBody['RightAileron']["CXtot"], avl.output["CXtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["CYtot"], avl.output["CYtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["CZtot"], avl.output["CZtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["Cltot"], avl.output["Cltot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["Cmtot"], avl.output["Cmtot"].deriv("RightAileron"), 8)
            self.assertAlmostEqual(ControlBody['RightAileron']["Cntot"], avl.output["Cntot"].deriv("RightAileron"), 8)

#==============================================================================
    def test_wing_tail(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 2.0

        avl.input.RollRate  = 0.1
        avl.input.PitchRate = 0.2
        avl.input.YawRate   = 0.3

        avl.input.CDp = 0.001

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        htail = {"groupName"         : "hTail",
                 "numChord"          : 8,
                 "spaceChord"        : 0.7,
                 "numSpanTotal"      : 10}

        vtail = {"numChord"          : 8,
                 "spaceChord"        : 1.0,
                 "numSpanTotal"      : 8,
                 "spaceSpan"         : 0.5}

        avl.input.AVL_Surface = {"Wing" : wing,
                                 "hTail": htail,
                                 "vTail": vtail}

        # print()
        # print("AlphaTrue  =", avl.output["Alpha" ].value)
        # print("BetaTrue   =", avl.output["Beta"  ].value)
        # print("MachTrue   =", avl.output["Mach"  ].value)
        # print("pb_2VTrue  =", avl.output["pb/2V" ].value)
        # print("qc_2VTrue  =", avl.output["qc/2V" ].value)
        # print("rb_2VTrue  =", avl.output["rb/2V" ].value)
        # print("pPb_2VTrue =", avl.output["p'b/2V"].value)
        # print("rPb_2VTrue =", avl.output["r'b/2V"].value)
        # print("CXtotTrue  =", avl.output["CXtot" ].value)
        # print("CYtotTrue  =", avl.output["CYtot" ].value)
        # print("CZtotTrue  =", avl.output["CZtot" ].value)
        # print("CltotTrue  =", avl.output["Cltot" ].value)
        # print("CmtotTrue  =", avl.output["Cmtot" ].value)
        # print("CntotTrue  =", avl.output["Cntot" ].value)
        # print("ClPtotTrue =", avl.output["Cl'tot"].value)
        # print("CnPtotTrue =", avl.output["Cn'tot"].value)
        # print("CLtotTrue  =", avl.output["CLtot" ].value)
        # print("CDtotTrue  =", avl.output["CDtot" ].value)
        # print("CDvisTrue  =", avl.output["CDvis" ].value)
        # print("CLffTrue   =", avl.output["CLff"  ].value)
        # print("CYffTrue   =", avl.output["CYff"  ].value)
        # print("CDffTrue   =", avl.output["CDff"  ].value)
        # print("CDindTrue  =", avl.output["CDind" ].value)
        # print("eTrue      =", avl.output["e"     ].value)

        AlphaTrue  = 1.0
        BetaTrue   = 2.0
        MachTrue   = 0.5
        pb_2VTrue  = 0.09474904758445407
        qc_2VTrue  = 0.2
        rb_2VTrue  = 0.3016995491906457
        pPb_2VTrue = 0.1
        rPb_2VTrue = 0.3
        CXtotTrue  = 0.2819749866379456
        CYtotTrue  = 0.3299042737982653
        CZtotTrue  = -1.859355248703562
        CltotTrue  = 0.07173904373817362
        CmtotTrue  = -2.008623540362779
        CntotTrue  = -0.1883105098861406
        ClPtotTrue = 0.06844164597939142
        CnPtotTrue = -0.1895338482321232
        CLtotTrue  = 1.863993201965148
        CDtotTrue  = -0.2494818169700332
        CDvisTrue  = 0.001
        CLffTrue   = 1.855347458479792
        CYffTrue   = 0.2529478491232741
        CDffTrue   = 0.6150573936973005
        CDindTrue  = -0.2504818169700332
        eTrue      = 0.604869804943736

        self.assertAlmostEqual(AlphaTrue , avl.output["Alpha" ].value, 4)
        self.assertAlmostEqual(BetaTrue  , avl.output["Beta"  ].value, 4)
        self.assertAlmostEqual(MachTrue  , avl.output["Mach"  ].value, 4)
        self.assertAlmostEqual(pb_2VTrue , avl.output["pb/2V" ].value, 4)
        self.assertAlmostEqual(qc_2VTrue , avl.output["qc/2V" ].value, 4)
        self.assertAlmostEqual(rb_2VTrue , avl.output["rb/2V" ].value, 4)
        self.assertAlmostEqual(pPb_2VTrue, avl.output["p'b/2V"].value, 4)
        self.assertAlmostEqual(rPb_2VTrue, avl.output["r'b/2V"].value, 4)
        self.assertAlmostEqual(CXtotTrue , avl.output["CXtot" ].value, 4)
        self.assertAlmostEqual(CYtotTrue , avl.output["CYtot" ].value, 4)
        self.assertAlmostEqual(CZtotTrue , avl.output["CZtot" ].value, 4)
        self.assertAlmostEqual(CltotTrue , avl.output["Cltot" ].value, 4)
        self.assertAlmostEqual(CmtotTrue , avl.output["Cmtot" ].value, 4)
        self.assertAlmostEqual(CntotTrue , avl.output["Cntot" ].value, 4)
        self.assertAlmostEqual(ClPtotTrue, avl.output["Cl'tot"].value, 4)
        self.assertAlmostEqual(CnPtotTrue, avl.output["Cn'tot"].value, 4)
        self.assertAlmostEqual(CLtotTrue , avl.output["CLtot" ].value, 4)
        self.assertAlmostEqual(CDtotTrue , avl.output["CDtot" ].value, 4)
        self.assertAlmostEqual(CDvisTrue , avl.output["CDvis" ].value, 4)
        self.assertAlmostEqual(CLffTrue  , avl.output["CLff"  ].value, 4)
        self.assertAlmostEqual(CYffTrue  , avl.output["CYff"  ].value, 4)
        self.assertAlmostEqual(CDffTrue  , avl.output["CDff"  ].value, 4)
        self.assertAlmostEqual(CDindTrue , avl.output["CDind" ].value, 4)
        self.assertAlmostEqual(eTrue     , avl.output["e"     ].value, 4)

        # Alpha stability derivatives
        # print("CLaTrue  =", avl.output["CLa" ].value)
        # print("CYaTrue  =", avl.output["CYa" ].value)
        # print("ClPaTrue =", avl.output["Cl'a"].value)
        # print("CmaTrue  =", avl.output["Cma" ].value)
        # print("CnPaTrue =", avl.output["Cn'a"].value)

        CLaTrue  = 4.118327586544591
        CYaTrue  = 0.3458124255025434
        ClPaTrue = -0.1626020019146856
        CmaTrue  = -0.9227102440442653
        CnPaTrue = -0.2071250653726279

        self.assertAlmostEqual(CLaTrue , avl.output["CLa" ].value, 4)
        self.assertAlmostEqual(CYaTrue , avl.output["CYa" ].value, 4)
        self.assertAlmostEqual(ClPaTrue, avl.output["Cl'a"].value, 4)
        self.assertAlmostEqual(CmaTrue , avl.output["Cma" ].value, 4)
        self.assertAlmostEqual(CnPaTrue, avl.output["Cn'a"].value, 4)

        # Beta stability derivatives
        # print("CLbTrue  =", avl.output["CLb" ].value)
        # print("CYbTrue  =", avl.output["CYb" ].value)
        # print("ClPbTrue =", avl.output["Cl'b"].value)
        # print("CmbTrue  =", avl.output["Cmb" ].value)
        # print("CnPbTrue =", avl.output["Cn'b"].value)

        CLbTrue  = -0.03655550147336076
        CYbTrue  = -0.4959910181047151
        ClPbTrue = -0.3207945915531573
        CmbTrue  = 0.06155514345290854
        CnPbTrue = 0.2992239603044715

        self.assertAlmostEqual(CLbTrue , avl.output["CLb" ].value, 4)
        self.assertAlmostEqual(CYbTrue , avl.output["CYb" ].value, 4)
        self.assertAlmostEqual(ClPbTrue, avl.output["Cl'b"].value, 4)
        self.assertAlmostEqual(CmbTrue , avl.output["Cmb" ].value, 4)
        self.assertAlmostEqual(CnPbTrue, avl.output["Cn'b"].value, 4)

        # Roll rate p' stability derivatives
        # print("CLpPTrue  =", avl.output["CLp'" ].value)
        # print("CYpPTrue  =", avl.output["CYp'" ].value)
        # print("ClPpPTrue =", avl.output["Cl'p'"].value)
        # print("CmpPTrue  =", avl.output["Cmp'" ].value)
        # print("CnPpPTrue =", avl.output["Cn'p'"].value)

        CLpPTrue  = -0.2938642220505157
        CYpPTrue  = 0.7098931332923271
        ClPpPTrue = -0.2687129813055149
        CmpPTrue  = 0.2946143866496476
        CnPpPTrue = -0.3704619477289713

        self.assertAlmostEqual(CLpPTrue , avl.output["CLp'" ].value, 4)
        self.assertAlmostEqual(CYpPTrue , avl.output["CYp'" ].value, 4)
        self.assertAlmostEqual(ClPpPTrue, avl.output["Cl'p'"].value, 4)
        self.assertAlmostEqual(CmpPTrue , avl.output["Cmp'" ].value, 4)
        self.assertAlmostEqual(CnPpPTrue, avl.output["Cn'p'"].value, 4)

        # Pitch rate q' stability derivatives
        # print("CLqPTrue  =", avl.output["CLq'" ].value)
        # print("CYqPTrue  =", avl.output["CYq'" ].value)
        # print("ClPqPTrue =", avl.output["Cl'q'"].value)
        # print("CmqPTrue  =", avl.output["Cmq'" ].value)
        # print("CnPqPTrue =", avl.output["Cn'q'"].value)

        CLqPTrue  = 9.250119191858278
        CYqPTrue  = 0.377116595595087
        ClPqPTrue = 0.4029373196844015
        CmqPTrue  = -11.7657123289711
        CnPqPTrue = -0.1349303788865462

        self.assertAlmostEqual(CLqPTrue , avl.output["CLq'" ].value, 4)
        self.assertAlmostEqual(CYqPTrue , avl.output["CYq'" ].value, 4)
        self.assertAlmostEqual(ClPqPTrue, avl.output["Cl'q'"].value, 4)
        self.assertAlmostEqual(CmqPTrue , avl.output["Cmq'" ].value, 4)
        self.assertAlmostEqual(CnPqPTrue, avl.output["Cn'q'"].value, 4)

        # Yaw rate r' stability derivatives
        # print("CLrPTrue  =", avl.output["CLr'" ].value)
        # print("CYrPTrue  =", avl.output["CYr'" ].value)
        # print("ClPrPTrue =", avl.output["Cl'r'"].value)
        # print("CmrPTrue  =", avl.output["Cmr'" ].value)
        # print("CnPrPTrue =", avl.output["Cn'r'"].value)

        CLrPTrue  = -0.2333570033839645
        CYrPTrue  = 0.8805003102526999
        ClPrPTrue = 0.4712656557162264
        CmrPTrue  = -0.1277542701416303
        CnPrPTrue = -0.5771230036264411

        self.assertAlmostEqual(CLrPTrue , avl.output["CLr'" ].value, 4)
        self.assertAlmostEqual(CYrPTrue , avl.output["CYr'" ].value, 4)
        self.assertAlmostEqual(ClPrPTrue, avl.output["Cl'r'"].value, 4)
        self.assertAlmostEqual(CmrPTrue , avl.output["Cmr'" ].value, 4)
        self.assertAlmostEqual(CnPrPTrue, avl.output["Cn'r'"].value, 4)

        # Axial vel (body axis)  stability derivatives
        # print("CXuTrue =", avl.output["CXu"].value)
        # print("CYuTrue =", avl.output["CYu"].value)
        # print("CZuTrue =", avl.output["CZu"].value)
        # print("CluTrue =", avl.output["Clu"].value)
        # print("CmuTrue =", avl.output["Cmu"].value)
        # print("CnuTrue =", avl.output["Cnu"].value)

        CXuTrue = -0.07028198570508082
        CYuTrue = 0.2625057597983784
        CZuTrue = -1.909949786271903
        CluTrue = -0.0470370424702159
        CmuTrue = -1.638247526986853
        CnuTrue = -0.1515868083855788
        
        self.assertAlmostEqual(CXuTrue, avl.output["CXu"].value, 4)
        self.assertAlmostEqual(CYuTrue, avl.output["CYu"].value, 4)
        self.assertAlmostEqual(CZuTrue, avl.output["CZu"].value, 4)
        self.assertAlmostEqual(CluTrue, avl.output["Clu"].value, 4)
        self.assertAlmostEqual(CmuTrue, avl.output["Cmu"].value, 4)
        self.assertAlmostEqual(CnuTrue, avl.output["Cnu"].value, 4)

        # Sideslip vel (body axis) stability derivatives
        # print("CXvTrue =", avl.output["CXv"].value)
        # print("CYvTrue =", avl.output["CYv"].value)
        # print("CZvTrue =", avl.output["CZv"].value)
        # print("ClvTrue =", avl.output["Clv"].value)
        # print("CmvTrue =", avl.output["Cmv"].value)
        # print("CnvTrue =", avl.output["Cnv"].value)

        CXvTrue = -0.1891232938404756
        CYvTrue = -0.4869903200321842
        CZvTrue = -0.03578160575317908
        ClvTrue = -0.3277141927214199
        CmvTrue = 0.00375061073859102
        CnvTrue = 0.2884140196188684

        self.assertAlmostEqual(CXvTrue, avl.output["CXv"].value, 4)
        self.assertAlmostEqual(CYvTrue, avl.output["CYv"].value, 4)
        self.assertAlmostEqual(CZvTrue, avl.output["CZv"].value, 4)
        self.assertAlmostEqual(ClvTrue, avl.output["Clv"].value, 4)
        self.assertAlmostEqual(CmvTrue, avl.output["Cmv"].value, 4)
        self.assertAlmostEqual(CnvTrue, avl.output["Cnv"].value, 4)

        # Normal vel (body axis)  stability derivatives
        # print("CXwTrue =", avl.output["CXw"].value)
        # print("CYwTrue =", avl.output["CYw"].value)
        # print("CZwTrue =", avl.output["CZw"].value)
        # print("ClwTrue =", avl.output["Clw"].value)
        # print("CmwTrue =", avl.output["Cmw"].value)
        # print("CnwTrue =", avl.output["Cnw"].value)

        CXwTrue = 0.6825521730032321
        CYwTrue = 0.225644885554448
        CZwTrue = -3.958668601535436
        ClwTrue = 0.1554346272123236
        CmwTrue = -1.053245857298092
        CnwTrue = -0.08525360021500156

        self.assertAlmostEqual(CXwTrue, avl.output["CXw"].value, 4)
        self.assertAlmostEqual(CYwTrue, avl.output["CYw"].value, 4)
        self.assertAlmostEqual(CZwTrue, avl.output["CZw"].value, 4)
        self.assertAlmostEqual(ClwTrue, avl.output["Clw"].value, 4)
        self.assertAlmostEqual(CmwTrue, avl.output["Cmw"].value, 4)
        self.assertAlmostEqual(CnwTrue, avl.output["Cnw"].value, 4)

        # Roll rate (body axis)  stability derivatives
        # print("CXpTrue =", avl.output["CXp"].value)
        # print("CYpTrue =", avl.output["CYp"].value)
        # print("CZpTrue =", avl.output["CZp"].value)
        # print("ClpTrue =", avl.output["Clp"].value)
        # print("CmpTrue =", avl.output["Cmp"].value)
        # print("CnpTrue =", avl.output["Cnp"].value)

        CXpTrue = 0.02734084195132762
        CYpTrue = 0.6944181638469976
        CZpTrue = 0.290268196576447
        ClpTrue = -0.2705659181690844
        CmpTrue = 0.2967991348981743
        CnpTrue = -0.365110973898108

        self.assertAlmostEqual(CXpTrue, avl.output["CXp"].value, 4)
        self.assertAlmostEqual(CYpTrue, avl.output["CYp"].value, 4)
        self.assertAlmostEqual(CZpTrue, avl.output["CZp"].value, 4)
        self.assertAlmostEqual(ClpTrue, avl.output["Clp"].value, 4)
        self.assertAlmostEqual(CmpTrue, avl.output["Cmp"].value, 4)
        self.assertAlmostEqual(CnpTrue, avl.output["Cnp"].value, 4)

        # Pitch rate (body axis)  stability derivatives
        # print("CXqTrue =", avl.output["CXq"].value)
        # print("CYqTrue =", avl.output["CYq"].value)
        # print("CZqTrue =", avl.output["CZq"].value)
        # print("ClqTrue =", avl.output["Clq"].value)
        # print("CmqTrue =", avl.output["Cmq"].value)
        # print("CnqTrue =", avl.output["Cnq"].value)

        CXqTrue = 2.720908021350824
        CYqTrue = 0.377116595595087
        CZqTrue = -9.204034618244297
        ClqTrue = 0.4052308101920075
        CmqTrue = -11.7657123289711
        CnqTrue = -0.12787760246441

        self.assertAlmostEqual(CXqTrue, avl.output["CXq"].value, 4)
        self.assertAlmostEqual(CYqTrue, avl.output["CYq"].value, 4)
        self.assertAlmostEqual(CZqTrue, avl.output["CZq"].value, 4)
        self.assertAlmostEqual(ClqTrue, avl.output["Clq"].value, 4)
        self.assertAlmostEqual(CmqTrue, avl.output["Cmq"].value, 4)
        self.assertAlmostEqual(CnqTrue, avl.output["Cnq"].value, 4)

        # Yaw rate (body axis)  stability derivatives
        # print("CXrTrue =", avl.output["CXr"].value)
        # print("CYrTrue =", avl.output["CYr"].value)
        # print("CZrTrue =", avl.output["CZr"].value)
        # print("ClrTrue =", avl.output["Clr"].value)
        # print("CmrTrue =", avl.output["Cmr"].value)
        # print("CnrTrue =", avl.output["Cnr"].value)

        CXrTrue = 0.2721347417070865
        CYrTrue = 0.8927555492799036
        CZrTrue = 0.2432365520456572
        ClrTrue = 0.4766166295470897
        CmrTrue = -0.1225930825294154
        CnrTrue = -0.5752700667628715

        self.assertAlmostEqual(CXrTrue, avl.output["CXr"].value, 4)
        self.assertAlmostEqual(CYrTrue, avl.output["CYr"].value, 4)
        self.assertAlmostEqual(CZrTrue, avl.output["CZr"].value, 4)
        self.assertAlmostEqual(ClrTrue, avl.output["Clr"].value, 4)
        self.assertAlmostEqual(CmrTrue, avl.output["Cmr"].value, 4)
        self.assertAlmostEqual(CnrTrue, avl.output["Cnr"].value, 4)

        # Neutral point and CG
        # print("XnpTrue =", avl.output["Xnp"].value)
        # print("XcgTrue =", avl.output["Xcg"].value)
        # print("YcgTrue =", avl.output["Ycg"].value)
        # print("ZcgTrue =", avl.output["Zcg"].value)

        XnpTrue = 0.4740497446242368
        XcgTrue = 0.25
        YcgTrue = 0.0
        ZcgTrue = 0.0

        self.assertAlmostEqual(XnpTrue, avl.output["Xnp"].value, 4)
        self.assertAlmostEqual(XcgTrue, avl.output["Xcg"].value, 4)
        self.assertAlmostEqual(YcgTrue, avl.output["Ycg"].value, 4)
        self.assertAlmostEqual(ZcgTrue, avl.output["Zcg"].value, 4)

        # ControlStability
        #print("ControlStabilityTrue =", avl.output["ControlStability"].value)

        ControlStabilityTrue = {'LeftAileron': {'CLtot': 0.01236470863456, 'CYtot': -0.0005711861518126, "Cl'tot": 0.003245812016112, 'Cmtot': -0.003234883458759, "Cn'tot": 0.0007714186936987, 'CDff': 0.001269960311541, 'e': 0.005601470784052}, 
                                'RightAileron': {'CLtot': 0.005704488876958, 'CYtot': 0.0006969654416042, "Cl'tot": -0.0012240864846, 'Cmtot': -0.0009944240502486, "Cn'tot": -0.0001620437721934, 'CDff': 0.001909597542066, 'e': 0.002431894262731}, 
                                'Elevator': {'CLtot': 0.009188221421924, 'CYtot': 0.0003305059144593, "Cl'tot": 0.0005072632713486, 'Cmtot': -0.02162289354974, "Cn'tot": -0.0003467959440848, 'CDff': 0.00779943700234, 'e': -0.00178014685464}, 
                                'Rudder': {'CLtot': -0.0002396650705303, 'CYtot': 0.004880169511806, "Cl'tot": -1.905382764638e-05, 'Cmtot': 0.0005653812602337, "Cn'tot": -0.003871221575955, 'CDff': 0.002048247085248, 'e': -0.001607133755587}}

        ControlStability = avl.output["ControlStability"].value

        for control in ControlStabilityTrue.keys():
            for coeff in ControlStabilityTrue[control].keys():
                self.assertAlmostEqual(ControlStabilityTrue[control][coeff], ControlStability[control][coeff], 4)

        # ControlBody
        #print("ControlBodyTrue =", avl.output["ControlBody"].value)

        ControlBodyTrue = {'LeftAileron': {'CXtot': 0.0005327280438995, 'CYtot': -0.0005711861518126, 'CZtot': -0.01235729332384, 'Cltot': 0.003245812016112, 'Cmtot': -0.003234883458759, 'Cntot': 0.0007714186936987}, 
                           'RightAileron': {'CXtot': -0.0003924356091698, 'CYtot': 0.0006969654416042, 'CZtot': -0.005712207819628, 'Cltot': -0.0012240864846, 'Cmtot': -0.0009944240502486, 'Cntot': -0.0001620437721934}, 
                           'Elevator': {'CXtot': 0.001116872349839, 'CYtot': 0.0003305059144593, 'CZtot': -0.009170125966337, 'Cltot': 0.0005072632713486, 'Cmtot': -0.02162289354974, 'Cntot': -0.0003467959440848}, 
                           'Rudder': {'CXtot': 0.0003165733098541, 'CYtot': 0.004880169511806, 'CZtot': 0.0002452273859197, 'Cltot': -1.905382764638e-05, 'Cmtot': 0.0005653812602337, 'Cntot': -0.003871221575955}}

        ControlBody = avl.output["ControlBody"].value

        for control in ControlBodyTrue.keys():
            for coeff in ControlBodyTrue[control].keys():
                self.assertAlmostEqual(ControlBodyTrue[control][coeff], ControlBody[control][coeff], 4)

        # HingeMoment
        #print("HingeMomentTrue =", avl.output["HingeMoment"].value)

        HingeMomentTrue = {'LeftAileron': -0.003436585479755, 'RightAileron': -0.0022599126024, 'Elevator': -0.000153131475531, 'Rudder': -0.002011277198011}

        HingeMoment = avl.output["HingeMoment"].value

        for control in HingeMomentTrue.keys():
            self.assertAlmostEqual(HingeMomentTrue[control], HingeMoment[control], 4)

        # StripForces
        #print("StripForcesTrue =", avl.output["StripForces"].value)

        StripForcesTrue = {'Wing': {'Xle': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                    'Yle': [-1.495722430687, -1.461939766258, -1.396676670151, -1.304380714515, -1.191341716198, -1.06526309613, -0.9347369039141, -0.8086582838434, -0.6956192855209, -0.603323329876, -0.5380602337594, -0.5042775693187, -0.4957224306871, -0.4619397662578, -0.3966766701513, -0.3043807145147, -0.1913417161981, -0.06526309613043, 0.06526309608591, 0.1913417161566, 0.3043807144791, 0.396676670124, 0.4619397662406, 0.4957224306813, 0.5042775693129, 0.5380602337422, 0.6033233298487, 0.6956192854853, 0.8086582838019, 0.9347369038696, 1.065263096086, 1.191341716157, 1.304380714479, 1.396676670124, 1.461939766241, 1.495722430681], 
                                    'Zle': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                    'Chord': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
                                    'Area': [0.0170370868545, 0.04995021124954, 0.07945931129475, 0.1035533905882, 0.1205904774436, 0.1294095225469, 0.1294095225484, 0.1205904774481, 0.1035533905954, 0.07945931130405, 0.04995021126037, 0.01703708686612, 0.0170370868545, 0.04995021124954, 0.07945931129475, 0.1035533905882, 0.1205904774436, 0.1294095225469, 0.1294095225484, 0.1205904774481, 0.1035533905954, 0.07945931130405, 0.04995021126037, 0.01703708686612, 0.0170370868545, 0.04995021124954, 0.07945931129475, 0.1035533905882, 0.1205904774436, 0.1294095225469, 0.1294095225484, 0.1205904774481, 0.1035533905954, 0.07945931130405, 0.04995021126037, 0.01703708686612], 
                                    'c_cl': [0.05713745336977, 0.1722244276742, 0.2917293935614, 0.4207102400565, 0.5597708600895, 0.7034357527641, 0.8418639659122, 0.9643019468149, 1.062392021241, 1.132225999036, 1.174651618374, 1.193807852727, 1.195975886428, 1.21224050585, 1.238087207979, 1.265464004143, 1.2873502396, 1.299229501085, 1.299748550667, 1.290605991278, 1.275272239459, 1.258233662491, 1.244079976114, 1.236288256997, 1.234183709534, 1.226408536983, 1.210164981068, 1.182703733554, 1.139396705142, 1.074813717774, 0.984028215045, 0.8638120892249, 0.7134007810971, 0.5346967440865, 0.3321498283106, 0.1127876089501], 
                                    'ai': [-0.1102662053317, -0.1170953499571, -0.1237296078527, -0.1240778916249, -0.1188977592912, -0.1062066027969, -0.06739369309278, 0.0454654409955, 0.2426124465743, 0.4330952836578, 0.5676261159348, 0.6366874073991, 0.652520285533, 0.7066899824471, 0.7872497983733, 0.8747889681352, 1.053753533426, 1.264122101246, 0.3296120002264, 0.401755837549, 0.3337896478421, 0.3360163271176, 0.3229427925372, 0.3004886495189, 0.2943823349309, 0.2709943871446, 0.2262601110551, 0.1685614281867, 0.1220423852342, 0.1352571650226, 0.1937507121215, 0.2620967583214, 0.3254518390271, 0.3799730171961, 0.4236203203982, 0.4490091927713], 
                                    'cl_norm': [0.03392298204356, 0.1033469553206, 0.178607247847, 0.2649856392231, 0.3653156419972, 0.4781336562615, 0.5975247151573, 0.7143920106539, 0.8185936283862, 0.9013766259823, 0.9573182083882, 0.9849033621415, 0.9897817454934, 1.015835388739, 1.063011738117, 1.125016191906, 1.195235736851, 1.267373556555, 1.335994639893, 1.397015503576, 1.447411799244, 1.485513962866, 1.510962388105, 1.523866661829, 1.526982190407, 1.539998699324, 1.564131712333, 1.593367268712, 1.616651681151, 1.617936646119, 1.577256557034, 1.473396198355, 1.288297820721, 1.012576542828, 0.6508177836297, 0.2249650029039], 
                                    'cl': [0.05754234631622, 0.1734487717268, 0.2936173642937, 0.4229017520349, 0.5619299662807, 0.7053165970252, 0.843316925323, 0.9652634605514, 1.062870602056, 1.132288360176, 1.174409217252, 1.193404048673, 1.195552921643, 1.21184009367, 1.237731190619, 1.265172313767, 1.287138775406, 1.299107870234, 1.299718644321, 1.290662471485, 1.275404199385, 1.258426760773, 1.244317725152, 1.23655096619, 1.234436091254, 1.226566135991, 1.210147981486, 1.182453125513, 1.138878101988, 1.074022735807, 0.9829976858008, 0.8626198942447, 0.712173914264, 0.5336041495042, 0.3313782844071, 0.1125060293309], 
                                    'cd': [-0.02269767036534, -0.06863964398363, -0.1056158869628, -0.1218801545092, -0.1188100645628, -0.1016147047376, -0.07589315336933, -0.04666974947172, -0.01814652177418, 0.006308119254701, 0.02413818035599, 0.03355212208636, 0.03466875970831, 0.03351866106348, 0.03120084425924, 0.02775446753931, 0.02334931091297, 0.01830643087219, 0.01305606714144, 0.008027199006105, 0.003569153124797, -8.215434605694e-05, -0.002763711793705, -0.004261675444456, -0.003688380992873, 0.001673864018737, 0.01153485521178, 0.02467862421227, 0.03965411925749, 0.05469508376063, 0.06762646498933, 0.075839178572, 0.07651290809484, 0.0672609294894, 0.04710034826269, 0.0171159681361], 
                                    'cdv': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                    'cm_c/4': [-0.06748884595161, -0.1892153641068, -0.279418754619, -0.335082573569, -0.3618214324275, -0.3685007693047, -0.3635937442864, -0.3535347912871, -0.3425544212373, -0.333161378267, -0.3266687568837, -0.3235030864916, -0.3228430574889, -0.3200176335216, -0.3149290468017, -0.3079736005342, -0.2993925162319, -0.2894509220469, -0.2786107566535, -0.2675788642496, -0.2571328664744, -0.2480950447953, -0.2413215461037, -0.2376228055137, -0.2368487085089, -0.2329418909414, -0.2250706891599, -0.213413390446, -0.1983963022365, -0.1805069573125, -0.1600935934053, -0.1372567250711, -0.1118504446629, -0.08357269975225, -0.05216032572956, -0.01780989809422], 
                                    'cm_LE': [-0.08177320929405, -0.2322714710254, -0.3523511030093, -0.4402601335832, -0.5017641474499, -0.5443597074957, -0.5740597357644, -0.5946102779908, -0.6081524265477, -0.616217878026, -0.6203316614771, -0.6219550496734, -0.6218370290958, -0.6230777599841, -0.6244508487964, -0.62433960157, -0.621230076132, -0.6142582973181, -0.6035478943202, -0.5902303620692, -0.575950926339, -0.5626534604181, -0.5523415401322, -0.546694869763, -0.5453946358924, -0.5395440251871, -0.5276119344269, -0.5090893238347, -0.4832454785221, -0.449210386756, -0.4061006471666, -0.3532097473773, -0.2902006399371, -0.2172468857739, -0.1351977828072, -0.04600680033175], 
                                    'C.P.x/c': [1.422855301741, 1.340900570947, 1.201642472819, 1.04234141726, 0.8938906165164, 0.7724615029036, 0.6811472156771, 0.6162573025246, 0.5722917451803, 0.5442372190553, 0.5281558183341, 0.5210759083241, 0.5200366095422, 0.5140757928321, 0.504440583859, 0.4934242333499, 0.4826031364703, 0.4728074578555, 0.4643623605546, 0.4573190087737, 0.4516089225662, 0.447146987436, 0.4439388479532, 0.4421657998827, 0.4418679388807, 0.4399138449254, 0.4359860881506, 0.4304835945218, 0.4242032811854, 0.4180662348148, 0.4128626350986, 0.4091161135824, 0.4070549586592, 0.4066192838453, 0.4074041757832, 0.4083017212512]}, 
                           'hTail': {'Xle': [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0], 
                                     'Yle': [-0.7316461936117, -0.5954194696178, -0.3750000000168, -0.1545805304094, -0.01835380639868, 0.01835380638828, 0.1545805303822, 0.3749999999832, 0.5954194695906, 0.7316461936013], 
                                     'Zle': [0.2926584774447, 0.2381677878471, 0.1500000000067, 0.06183221216376, 0.007341522559473, 0.007341522555311, 0.06183221215286, 0.1499999999933, 0.2381677878362, 0.2926584774405], 
                                     'Chord': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
                                     'Area': [0.0771356220447, 0.201943680258, 0.2496161164303, 0.2019436802712, 0.07713562206601, 0.0771356220447, 0.201943680258, 0.2496161164303, 0.2019436802712, 0.07713562206601], 
                                     'c_cl': [0.6224061229963, 1.62379884856, 2.135704125158, 2.369322576072, 3.036720862448, 2.980090250539, 2.073232328621, 1.539036528982, 0.9213001567941, 0.3059173542654], 
                                     'ai': [1.133745819631, 1.091276840926, 1.061854116598, 1.134387322159, 1.361096509151, 0.2294871928161, 0.4328666361631, 0.5283080495509, 0.5638231257371, 0.588791707622], 
                                     'cl_norm': [0.3516148588636, 0.9164138503373, 1.201345833181, 1.327109882014, 1.693994860963, 2.174822992017, 1.61291555908, 1.326314893891, 0.8784163973436, 0.3099741614128], 
                                     'cl': [0.6236043216807, 1.629189411825, 2.144584757115, 2.379680566235, 3.046425903536, 2.987889642148, 2.081063911847, 1.545940882195, 0.9258990224991, 0.3073960148784], 
                                     'cd': [-0.06909264939192, -0.3199746166153, -0.5313187622793, -0.6206111026698, -0.5748937805069, -0.4564322095758, -0.4658224979551, -0.4130414388361, -0.2759547592906, -0.08863427262797], 
                                     'cdv': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                     'cm_c/4': [0.047567434489, 0.09217685574935, 0.08399206984827, 0.08125644335263, 0.177054500232, 0.175376680813, 0.06064876053577, 0.03316353420203, 0.02419325892098, 0.0115530219463],
                                     'cm_LE': [-0.1080340962601, -0.3137728563907, -0.4499339614413, -0.5110742006654, -0.58212571538, -0.5696458818218, -0.4576593216196, -0.3515955980434, -0.2061317802776, -0.06492631662004], 
                                     'C.P.x/c': [0.1737217690204, 0.193421645709, 0.2108352761206, 0.2158540543191, 0.1918812385929, 0.191304163869, 0.2208568486578, 0.2285479932745, 0.223870521155, 0.2124164875694]}, 
                           'vTail': {'Xle': [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0], 
                                     'Yle': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                     'Zle': [0.9903926402022, 0.9157348061559, 0.7777851165215, 0.5975451610273, 0.4024548390167, 0.2222148835159, 0.08426519386899, 0.009607359806596], 
                                     'Chord': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
                                     'Area': [0.03806023374221, 0.1083863756566, 0.1622116744031, 0.1913417161757, 0.1913417161791, 0.1622116744128, 0.1083863756711, 0.03806023375939], 
                                     'c_cl': [0.1480011067087, 0.4317707690465, 0.672558462148, 0.8225401096563, 0.8262184774058, 0.6555259285669, 0.3748337127576, 0.1153618532726], 
                                     'ai': [0.6250615325411, 0.5930071486616, 0.540862233547, 0.47295305416, 0.3981906486921, 0.3429181207828, 0.3346642386895, 0.3690300254941], 
                                     'cl_norm': [0.3243056612754, 0.870211170169, 1.170629503485, 1.198490308722, 1.00842211554, 0.6876443403844, 0.3526089134319, 0.1025418567408], 
                                     'cl': [0.1480011067087, 0.4317707690465, 0.672558462148, 0.8225401096563, 0.8262184774058, 0.6555259285669, 0.3748337127576, 0.1153618532726], 
                                     'cd': [-0.0007312986365578, -0.01876079100566, -0.05656985315704, -0.09920081652926, -0.1246300396941, -0.1155978618405, -0.07387961755708, -0.02237394311633], 
                                     'cdv': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                                     'cm_c/4': [-0.002423459151371, -0.006206411349657, -0.005759923979807, 0.002480876742362, 0.01791748285507, 0.03076787929747, 0.0283260185829, 0.01161136237024], 
                                     'cm_LE': [-0.03942373582855, -0.1141491036113, -0.1738995395168, -0.2031541506717, -0.1886371364964, -0.1331136028442, -0.0653824096065, -0.0172291009479], 
                                     'C.P.x/c': [0.2663746015504, 0.2643743203445, 0.2585641982132, 0.246983883566, 0.2283138681293, 0.2030638256144, 0.1744304404358, 0.1493483370728]}
                           }

        StripForces = avl.output["StripForces"].value

        for surf in StripForcesTrue.keys():
            for name in StripForcesTrue[surf].keys():
                for i in range(len(StripForcesTrue[surf][name])):
                    self.assertAlmostEqual(StripForcesTrue[surf][name][i], StripForces[surf][name][i], 3)

#==============================================================================
    def test_wing_Vtail(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        vtail = {"numChord"          : 8,
                 "spaceChord"        : 1.0,
                 "numSpanTotal"      : 8,
                 "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing"  : wing,
                                 "VTail1": vtail,
                                 "VTail2": vtail}

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot
        #print("CLtot = ", CLtot)
        #print("CDtot = ", CDtot)

        CLtotTrue  = 0.14681
        CDtotTrue  = 0.02022
        self.assertAlmostEqual(CLtotTrue, CLtot, 4)
        self.assertAlmostEqual(CDtotTrue, CDtot, 4)


#==============================================================================
    def test_wing_tail_control(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 2.0

        avl.input.RollRate  = 0.1
        avl.input.PitchRate = 0.2
        avl.input.YawRate   = 0.3

        avl.input.CDp = 0.001

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        htail = {"groupName"         : "hTail",
                 "numChord"          : 8,
                 "spaceChord"        : 0.7,
                 "numSpanTotal"      : 10}

        vtail = {"numChord"          : 8,
                 "spaceChord"        : 1.0,
                 "numSpanTotal"      : 8,
                 "spaceSpan"         : 0.5}

        avl.input.AVL_Surface = {"Wing" : wing,
                                 "hTail": htail,
                                 "vTail": vtail}

        # Set control surface parameters
        aileronLeft  = {"deflectionAngle" : -25.0}
        aileronRight = {"deflectionAngle" :  25.0}
        elevator     = {"deflectionAngle" :  5.0}
        rudder       = {"deflectionAngle" : -2.0}

        avl.input.AVL_Control = {"LeftAileron" : aileronLeft ,
                                 "RightAileron": aileronRight,
                                 "Elevator"    : elevator    ,
                                 "Rudder"      : rudder      }


        # print()
        # print("AlphaTrue  =", avl.output["Alpha" ].value)
        # print("BetaTrue   =", avl.output["Beta"  ].value)
        # print("MachTrue   =", avl.output["Mach"  ].value)
        # print("pb_2VTrue  =", avl.output["pb/2V" ].value)
        # print("qc_2VTrue  =", avl.output["qc/2V" ].value)
        # print("rb_2VTrue  =", avl.output["rb/2V" ].value)
        # print("pPb_2VTrue =", avl.output["p'b/2V"].value)
        # print("rPb_2VTrue =", avl.output["r'b/2V"].value)
        # print("CXtotTrue  =", avl.output["CXtot" ].value)
        # print("CYtotTrue  =", avl.output["CYtot" ].value)
        # print("CZtotTrue  =", avl.output["CZtot" ].value)
        # print("CltotTrue  =", avl.output["Cltot" ].value)
        # print("CmtotTrue  =", avl.output["Cmtot" ].value)
        # print("CntotTrue  =", avl.output["Cntot" ].value)
        # print("ClPtotTrue =", avl.output["Cl'tot"].value)
        # print("CnPtotTrue =", avl.output["Cn'tot"].value)
        # print("CLtotTrue  =", avl.output["CLtot" ].value)
        # print("CDtotTrue  =", avl.output["CDtot" ].value)
        # print("CDvisTrue  =", avl.output["CDvis" ].value)
        # print("CLffTrue   =", avl.output["CLff"  ].value)
        # print("CYffTrue   =", avl.output["CYff"  ].value)
        # print("CDffTrue   =", avl.output["CDff"  ].value)
        # print("CDindTrue  =", avl.output["CDind" ].value)
        # print("eTrue      =", avl.output["e"     ].value)

        AlphaTrue  = 1.0
        BetaTrue   = 2.0
        MachTrue   = 0.5
        pb_2VTrue  = 0.09474904758445407
        qc_2VTrue  = 0.2
        rb_2VTrue  = 0.3016995491906457
        pPb_2VTrue = 0.1
        rPb_2VTrue = 0.3
        CXtotTrue  = 0.2434511334638081
        CYtotTrue  = 0.3534396878923751
        CZtotTrue  = -1.74137092030149
        CltotTrue  = -0.03771690397564934
        CmtotTrue  = -2.069696065509661
        CntotTrue  = -0.1986037048209423
        ClPtotTrue = -0.04117727208497329
        CnPtotTrue = -0.1979152057770004
        CLtotTrue  = 1.745354509204636
        CDtotTrue  = -0.2130229416178314
        CDvisTrue  = 0.001
        CLffTrue   = 1.797634510552186
        CYffTrue   = 0.2766279456798953
        CDffTrue   = 0.705579146387139
        CDindTrue  = -0.2140229416178314
        eTrue      = 0.4974510184993075

        self.assertAlmostEqual(AlphaTrue , avl.output["Alpha" ].value, 4)
        self.assertAlmostEqual(BetaTrue  , avl.output["Beta"  ].value, 4)
        self.assertAlmostEqual(MachTrue  , avl.output["Mach"  ].value, 4)
        self.assertAlmostEqual(pb_2VTrue , avl.output["pb/2V" ].value, 4)
        self.assertAlmostEqual(qc_2VTrue , avl.output["qc/2V" ].value, 4)
        self.assertAlmostEqual(rb_2VTrue , avl.output["rb/2V" ].value, 4)
        self.assertAlmostEqual(pPb_2VTrue, avl.output["p'b/2V"].value, 4)
        self.assertAlmostEqual(rPb_2VTrue, avl.output["r'b/2V"].value, 4)
        self.assertAlmostEqual(CXtotTrue , avl.output["CXtot" ].value, 4)
        self.assertAlmostEqual(CYtotTrue , avl.output["CYtot" ].value, 4)
        self.assertAlmostEqual(CZtotTrue , avl.output["CZtot" ].value, 4)
        self.assertAlmostEqual(CltotTrue , avl.output["Cltot" ].value, 4)
        self.assertAlmostEqual(CmtotTrue , avl.output["Cmtot" ].value, 4)
        self.assertAlmostEqual(CntotTrue , avl.output["Cntot" ].value, 4)
        self.assertAlmostEqual(ClPtotTrue, avl.output["Cl'tot"].value, 4)
        self.assertAlmostEqual(CnPtotTrue, avl.output["Cn'tot"].value, 4)
        self.assertAlmostEqual(CLtotTrue , avl.output["CLtot" ].value, 4)
        self.assertAlmostEqual(CDtotTrue , avl.output["CDtot" ].value, 4)
        self.assertAlmostEqual(CDvisTrue , avl.output["CDvis" ].value, 4)
        self.assertAlmostEqual(CLffTrue  , avl.output["CLff"  ].value, 4)
        self.assertAlmostEqual(CYffTrue  , avl.output["CYff"  ].value, 4)
        self.assertAlmostEqual(CDffTrue  , avl.output["CDff"  ].value, 4)
        self.assertAlmostEqual(CDindTrue , avl.output["CDind" ].value, 4)
        self.assertAlmostEqual(eTrue     , avl.output["e"     ].value, 4)

#==============================================================================
    def test_geom_change(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing",
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot
        #print("CLtot = ", CLtot)
        #print("CDtot = ", CDtot)

        CLtotTrue  = 0.29683
        CDtotTrue  = 0.01617
        self.assertAlmostEqual(CLtotTrue, CLtot, 4)
        self.assertAlmostEqual(CDtotTrue, CDtot, 4)

        # change the geometry and rerun
        camber = self.myProblem.geometry.despmtr.camber
        self.myProblem.geometry.despmtr.camber = 1.1*camber

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot
        #print("CLtot = ", CLtot)
        #print("CDtot = ", CDtot)

        CLtotTrue  = 0.32434
        CDtotTrue  = 0.01801
        self.assertAlmostEqual(CLtotTrue, CLtot, 4)
        self.assertAlmostEqual(CDtotTrue, CDtot, 4)

        # reset the geometry
        self.myProblem.geometry.despmtr.camber = camber

#==============================================================================
    def test_MassProp_noUnits(self):

        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        mass = 0.1773
        x =  0.02463
        y = 0.
        z = 0.2239
        Ixx = 1.350
        Iyy = 0.7509
        Izz = 2.095

        avl.input.MassProp = {"Aircraft": {"mass":mass, "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz]}}

        # check there are errors if information is missing
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.Gravity  = 9.8         # m/s
        avl.input.Density  = 1.22557083  # kg/m^3
        avl.input.Velocity = 19.67167008 # "m/s"

        # make sure there are no errsos
        avl.runAnalysis()

        avl.input.MassProp = {"Aircraft":{"mass":mass, "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz, 1.0, 2.0, 3.0]},
                              "Engine"  :{"mass":mass, "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz]}}

        # again should not cause errors
        avl.runAnalysis()

        # test error handling of the mass properties parsing

        avl.input.MassProp = {"Aircraft": "1"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": "foo"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":mass, "CG":[x,y,z]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

        avl.input.MassProp = {"Aircraft": {"mass":"mass", "CG":[x,y,z], "massInertia":[Ixx, Iyy, Izz]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")


#==============================================================================
    def test_MassProp_Units(self):

        m    = pyCAPS.Unit("meter")
        kg   = pyCAPS.Unit("kg")
        s    = pyCAPS.Unit("s")
        K    = pyCAPS.Unit("Kelvin")
        deg  = pyCAPS.Unit("degree")
        ft   = pyCAPS.Unit("ft")
        slug = pyCAPS.Unit("slug")


        # Load avl aim
        avl = self.myProblem.analysis.create(aim = "avlAIM",
                                             unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0 * deg
        avl.input.Beta  = 0.0 * deg

        wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        mass = 0.1773
        x =  0.02463
        y = 0.
        z = 0.2239
        Ixx = 1.350
        Iyy = 0.7509
        Izz = 2.095

        avl.input.MassProp = {"Aircraft": {"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz] * kg*m**2}}

        # check there are errors if information is missing
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.Gravity  = 32.18 * ft/s**2
        avl.input.Density  = 0.002378 * slug/ft**3
        avl.input.Velocity = 64.5396 * ft/s

        # make sure there are no errsos
        avl.runAnalysis()

        avl.input.MassProp = {"Aircraft":{"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz, 1.0, 2.0, 3.0] * kg*m**2},
                              "Engine"  :{"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz] * kg*m**2}}

        # again should not cause errors
        avl.runAnalysis()

        # test error handling of the mass properties parsing

        avl.input.MassProp = {"Aircraft": "1"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": "foo"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":[mass,"kg"], "CG":[[x,y,z],"m"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

        avl.input.MassProp = {"Aircraft": {"mass":(mass), "CG":[[x,y,z],"m"], "massInertia":[[Ixx, Iyy, Izz], "kg*m^2"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":[mass,"kg"], "CG":[x,y,z], "massInertia":[[Ixx, Iyy, Izz], "kg*m^2"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":[mass,"kg"], "CG":[[x,y,z],"m"], "massInertia":[Ixx, Iyy, Izz]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        avl.input.MassProp = {"Aircraft": {"mass":("mass","kg"), "CG":[[x,y,z],"m"], "massInertia":[[Ixx, Iyy, Izz], "kg*m^2"]}}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            avl.runAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    def test_phase(self):

        # Initialize Problem object
        problemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile="../csmData/avlSections.csm", outLevel=0)

        # Load avl aim
        avl = myProblem.analysis.create(aim = "avlAIM", name = "avl")

        # Set new Mach/Alt parameters
        avl.input.Mach  = 0.5
        avl.input.Alpha = 1.0
        avl.input.Beta  = 0.0

        wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
                "numChord"          : 8,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}

        avl.input.AVL_Surface = {"Wing": wing}

        # Retrieve results
        CLtot = avl.output.CLtot
        CDtot = avl.output.CDtot

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0)

        avl = myProblem.analysis["avl"]

        # Retrieve results
        self.assertAlmostEqual(CLtot, avl.output.CLtot)
        self.assertAlmostEqual(CDtot, avl.output.CDtot)

        avl.input.Alpha  = 3.0
        self.assertAlmostEqual(0.41257, avl.output.CLtot, 4)


if __name__ == '__main__':
    unittest.main()
