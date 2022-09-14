import pyCAPS
import os
import numpy as np
import shutil
import math
import copy
from corsairlite.analysis.wrappers.CAPS.MSES.readSensx import readSensx

def msesRunNACA4_Alpha(camber,maxloc,thickness,Mach,Re,Alpha,
                         Acrit = None, GridAlpha = None, 
                         Coarse_Iteration = None, Fine_Iteration = None,
                         xGridRange = None, yGridRange = None, cleanup = True ):
    
    return msesRunNACA4( camber, maxloc, thickness,
                         Mach, Re,        
                         Alpha = Alpha,
                         CL = None,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )
    
def msesRunNACA4_CL(camber,maxloc,thickness,Mach,Re,CL,
                         Acrit = None, GridAlpha = None, 
                         Coarse_Iteration = 50, Fine_Iteration = None,
                         xGridRange = None, yGridRange = None, cleanup = True ):
    
    return msesRunNACA4( camber, maxloc, thickness,
                         Mach, Re,        
                         Alpha = None,
                         CL = CL,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )



def msesRunNACA4(camber, maxloc, thickness,
                 Mach, Re,        
                 Alpha = None,
                 CL = None,
                 Acrit = None, GridAlpha = None, 
                 Coarse_Iteration = None, Fine_Iteration = None,
                 xGridRange = None, yGridRange = None, cleanup = True ):

    geometryScript = "NACA4_airfoilSection.csm"
    runDirectory = "msesAnalysisTest"

    pstr = ""
    pstr += "attribute capsIntent $LINEARAERO \n"
    pstr += "attribute capsAIM $xfoilAIM;tsfoilAIM;msesAIM \n"
    pstr += "despmtr   camber       %f \n"%(camber)
    pstr += "despmtr   maxloc       %f \n"%(maxloc)
    pstr += "despmtr   thickness    %f \n"%(thickness)
    pstr += "udparg naca thickness  thickness \n"
    pstr += "udparg naca camber     camber \n"
    pstr += "udparg naca maxloc     maxloc \n"
    pstr += "udparg naca sharpte    1 \n"
    pstr += "udprim naca \n"
    pstr += "extract 0 \n"

    workDir = os.path.join(".",runDirectory)
    if os.path.exists(workDir) and os.path.isdir(workDir):
        shutil.rmtree(workDir)
    os.mkdir(workDir)
    geometryScript_path = os.path.join(workDir,geometryScript)
    f = open(geometryScript_path,'w')
    f.write(pstr)
    f.close()

    myProblem = pyCAPS.Problem(problemName=workDir, capsFile=geometryScript_path, outLevel=0)
    mses = myProblem.analysis.create(aim = "msesAIM")
    
    mses.input.Mach = Mach
    mses.input.Re   = Re
    
    if CL is not None and Alpha is not None:
        raise ValueError('Conflict Detected:  Alpha and CL have both been set')
        
    if CL is None and Alpha is None:
        raise ValueError('Must specifiy either Alpha or CL')
        
    if CL is not None:
        if CL >= 3.0:
            raise ValueError('CL exceeds 3.0, did you mean to set Alpha?')
        mses.input.CL = CL
        
    if Alpha is not None:
        mses.input.Alpha = Alpha
        
    if Acrit is not None:
        mses.input.Acrit            = Acrit # Default 9.0
    if GridAlpha is not None:
        mses.input.GridAlpha        = GridAlpha # Default 0.0
    elif Alpha is not None:
        mses.input.GridAlpha        = Alpha # Default 0.0
    else:
        mses.input.GridAlpha        = float(math.ceil(180.0 * CL / (4*np.pi**2)))  # Approximating alpha wiht thin airfoil theory, best thought I have
    if Coarse_Iteration is not None:
        mses.input.Coarse_Iteration = Coarse_Iteration # Default 20
    if Fine_Iteration is not None:
        mses.input.Fine_Iteration   = Fine_Iteration   # Default 200
    if xGridRange is not None:
        mses.input.xGridRange       = xGridRange       # Default [-1.75, 2.75]
    if yGridRange is not None:
        mses.input.yGridRange       = yGridRange       # Default [-2.0, 2.5]
        
    mses.input.Airfoil_Points = 151 # points

    if Mach <= 0.1 :
        mses.input.ISMOM = 2 # see mses documentation for M<=0.1
    
    mses.input.Design_Variable = {"camber": {}, "maxloc":{}, "thickness": {}} # should correspond to CSM despmtr values
    # mses.geometry.despmtr.camber = camber
    # mses.geometry.despmtr.maxloc = maxloc
    # mses.geometry.despmtr.thickness = thickness

    delta_camber    = 0.0005
    delta_maxloc    = 0.005
    delta_thickness = 0.005
    delta_mach      = 0.001
    delta_re        = 100
    delta_alpha     = 0.07
    delta_cl        = 0.01
    deltaDict = {}
    deltaDict['camber']    = delta_camber
    deltaDict['maxloc']    = delta_maxloc
    deltaDict['thickness'] = delta_thickness
    deltaDict['Mach']      = delta_mach
    deltaDict['Re']        = delta_re
    deltaDict['Alpha']     = delta_alpha
    deltaDict['CL']        = delta_cl
    centralDiff = [True]*6

    inputList = []
    if Alpha is not None:
        inputStack0 = [camber, maxloc, thickness, Mach, Re, Alpha]
        inputList.append(inputStack0)
    if CL is not None:
        inputStack0 = [camber, maxloc, thickness, Mach, Re, CL]
        inputList.append(inputStack0)
    inputStack = copy.deepcopy(inputStack0)

    # Camber
    camber_l = camber - delta_camber
    camber_u = camber + delta_camber
    if camber_l < 0.0:
        inputList.append(None)
        centralDiff[0] = False
    else:
        inputStack = copy.deepcopy(inputStack0)
        inputStack[0] = camber_l
        inputList.append(inputStack)
    inputStack = copy.deepcopy(inputStack0)
    inputStack[0] = camber_u
    inputList.append(inputStack)

    # Maxloc, assumes not near 0 or 1, otherwise would need some checks
    maxloc_l = maxloc - delta_maxloc
    maxloc_u = maxloc + delta_maxloc
    inputStack = copy.deepcopy(inputStack0)
    inputStack[1] = maxloc_l
    inputList.append(inputStack)
    inputStack = copy.deepcopy(inputStack0)
    inputStack[1] = maxloc_u
    inputList.append(inputStack)

    # Thickness, assumes not near zero, otherwise checks would be needed
    thickness_l = thickness - delta_thickness
    thickness_u = thickness + delta_thickness
    inputStack = copy.deepcopy(inputStack0)
    inputStack[2] = thickness_l
    inputList.append(inputStack)
    inputStack = copy.deepcopy(inputStack0)
    inputStack[2] = thickness_u
    inputList.append(inputStack)

    # Mach
    mach_l = Mach - delta_mach
    mach_u = Mach + delta_mach
    if mach_l < 0.0:
        inputList.append(None)
        centralDiff[3] = False
    else:
        inputStack = copy.deepcopy(inputStack0)
        inputStack[3] = mach_l
        inputList.append(inputStack)
    inputStack = copy.deepcopy(inputStack0)
    inputStack[3] = mach_u
    inputList.append(inputStack)

    # Re
    re_l = Re - delta_re
    re_u = Re + delta_re
    inputStack = copy.deepcopy(inputStack0)
    inputStack[4] = re_l
    inputList.append(inputStack)
    inputStack = copy.deepcopy(inputStack0)
    inputStack[4] = re_u
    inputList.append(inputStack)

    if Alpha is not None:
        alpha_l = Alpha - delta_alpha
        alpha_u = Alpha + delta_alpha
        inputStack = copy.deepcopy(inputStack0)
        inputStack[5] = alpha_l
        inputList.append(inputStack)
        inputStack = copy.deepcopy(inputStack0)
        inputStack[5] = alpha_u
        inputList.append(inputStack)
    if CL is not None:
        cl_l = CL - delta_cl
        cl_u = CL + delta_cl
        inputStack = copy.deepcopy(inputStack0)
        inputStack[5] = cl_l
        inputList.append(inputStack)
        inputStack = copy.deepcopy(inputStack0)
        inputStack[5] = cl_u
        inputList.append(inputStack)

    def getOutputs(inputValues):
        mses.geometry.despmtr.camber    = inputValues[0]
        mses.geometry.despmtr.maxloc    = inputValues[1]
        mses.geometry.despmtr.thickness = inputValues[2]
        mses.input.Mach              = inputValues[3]
        mses.input.Re                = inputValues[4]

        if Alpha is not None:
            mses.input.Alpha = inputValues[5]
        if CL is not None:
            mses.input.CL = inputValues[5]
        
        firstOutput = mses.output.CD # this runs mses

        msesOutput = os.path.join(".",runDirectory,"Scratch","msesAIM0","msesOutput.txt")
        f = open(msesOutput,'r')
        outputText = f.read()
        f.close()
        outputLines = outputText.split('\n')
        lastSection = '\n'.join(outputLines[-32:])
        if 'Converged on tolerance' not in lastSection:
            print('Case may not be converged')
            # print("Mach: %f CL: %f Re: %e camber: %f maxloc: %f thickness: %f "%(Mach,CL,Re,camber,maxloc,thickness))

        outputValues = {}
        outputDerivatives = {}

        if Alpha is not None:
            outputValues["CL"]                = mses.output.CL
        if CL is not None:
            outputValues["Alpha"]             = mses.output["Alpha"].value

        # mses.output.items()
        # mses.output.keys()
            
        outputValues["CD"]                = mses.output.CD
        outputValues["CD_p"]              = mses.output.CD_p
        outputValues["CD_v"]              = mses.output.CD_v
        outputValues["CD_w"]              = mses.output.CD_w
        outputValues["CM"]                = mses.output.CM

        outputDerivatives["CD"] = {}
        outputDerivatives["CD"]["Mach"]      = mses.output["CD"].deriv("Mach")
        outputDerivatives["CD"]["Re"]        = mses.output["CD"].deriv("Re")
        outputDerivatives["CD"]["camber"]    = mses.output["CD"].deriv("camber")
        outputDerivatives["CD"]["maxloc"]    = mses.output["CD"].deriv("maxloc")
        outputDerivatives["CD"]["thickness"] = mses.output["CD"].deriv("thickness")

        outputDerivatives["CD_p"] = {}
        outputDerivatives["CD_p"]["Mach"]      = mses.output["CD_p"].deriv("Mach")
        outputDerivatives["CD_p"]["Re"]        = mses.output["CD_p"].deriv("Re")
        outputDerivatives["CD_p"]["camber"]    = mses.output["CD_p"].deriv("camber")
        outputDerivatives["CD_p"]["maxloc"]    = mses.output["CD_p"].deriv("maxloc")
        outputDerivatives["CD_p"]["thickness"] = mses.output["CD_p"].deriv("thickness")

        outputDerivatives["CD_v"] = {}
        outputDerivatives["CD_v"]["Mach"]      = mses.output["CD_v"].deriv("Mach")
        outputDerivatives["CD_v"]["Re"]        = mses.output["CD_v"].deriv("Re")
        outputDerivatives["CD_v"]["camber"]    = mses.output["CD_v"].deriv("camber")
        outputDerivatives["CD_v"]["maxloc"]    = mses.output["CD_v"].deriv("maxloc")
        outputDerivatives["CD_v"]["thickness"] = mses.output["CD_v"].deriv("thickness")

        outputDerivatives["CD_w"] = {}
        outputDerivatives["CD_w"]["Mach"]      = mses.output["CD_w"].deriv("Mach")
        outputDerivatives["CD_w"]["Re"]        = mses.output["CD_w"].deriv("Re")
        outputDerivatives["CD_w"]["camber"]    = mses.output["CD_w"].deriv("camber")
        outputDerivatives["CD_w"]["maxloc"]    = mses.output["CD_w"].deriv("maxloc")
        outputDerivatives["CD_w"]["thickness"] = mses.output["CD_w"].deriv("thickness")

        outputDerivatives["CM"] = {}
        outputDerivatives["CM"]["Mach"]      = mses.output["CM"].deriv("Mach")
        outputDerivatives["CM"]["Re"]        = mses.output["CM"].deriv("Re")
        outputDerivatives["CM"]["camber"]    = mses.output["CM"].deriv("camber")
        outputDerivatives["CM"]["maxloc"]    = mses.output["CM"].deriv("maxloc")
        outputDerivatives["CM"]["thickness"] = mses.output["CM"].deriv("thickness")
        
        outputDerivatives["cheby"] = {}
        outputDerivatives["cheby"]["thickness"] = mses.output["Cheby_Modes"].deriv("thickness")
        outputDerivatives["cheby"]["camber"]    = mses.output["Cheby_Modes"].deriv("camber")
        outputDerivatives["cheby"]["maxloc"]    = mses.output["Cheby_Modes"].deriv("maxloc")

        if Alpha is not None:
            outputDerivatives["CL"] = {}
            outputDerivatives["CL"]["Alpha"]     = mses.output["CL"].deriv("Alpha")
            outputDerivatives["CL"]["Mach"]      = mses.output["CL"].deriv("Mach")
            outputDerivatives["CL"]["Re"]        = mses.output["CL"].deriv("Re")
            outputDerivatives["CL"]["camber"]    = mses.output["CL"].deriv("camber")
            outputDerivatives["CL"]["maxloc"]    = mses.output["CL"].deriv("maxloc")
            outputDerivatives["CL"]["thickness"] = mses.output["CL"].deriv("thickness")
            
            outputDerivatives["CD"]["Alpha"]     = mses.output["CD"].deriv("Alpha")
            outputDerivatives["CD_p"]["Alpha"]   = mses.output["CD_p"].deriv("Alpha")
            outputDerivatives["CD_v"]["Alpha"]   = mses.output["CD_v"].deriv("Alpha")
            outputDerivatives["CD_w"]["Alpha"]   = mses.output["CD_w"].deriv("Alpha")
            outputDerivatives["CM"]["Alpha"]     = mses.output["CM"].deriv("Alpha")
            
        if CL is not None:
            outputDerivatives["Alpha"] = {}
            outputDerivatives["Alpha"]["CL"]        = mses.output["Alpha"].deriv("CL")
            outputDerivatives["Alpha"]["Mach"]      = mses.output["Alpha"].deriv("Mach")
            outputDerivatives["Alpha"]["Re"]        = mses.output["Alpha"].deriv("Re")
            outputDerivatives["Alpha"]["camber"]    = mses.output["Alpha"].deriv("camber")
            outputDerivatives["Alpha"]["maxloc"]    = mses.output["Alpha"].deriv("maxloc")
            outputDerivatives["Alpha"]["thickness"] = mses.output["Alpha"].deriv("thickness")
            
            outputDerivatives["CD"]["CL"]     = mses.output["CD"].deriv("CL")
            outputDerivatives["CD_p"]["CL"]   = mses.output["CD_p"].deriv("CL")
            outputDerivatives["CD_v"]["CL"]   = mses.output["CD_v"].deriv("CL")
            outputDerivatives["CD_w"]["CL"]   = mses.output["CD_w"].deriv("CL")
            outputDerivatives["CM"]["CL"]     = mses.output["CM"].deriv("CL")

        return([outputValues, outputDerivatives])


    outputList = []
    for i in range(0,len(inputList)):
    # for inputValues in inputList:
        inputValues = inputList[i]
        # print(inputValues)
        if inputValues is not None:
            try:
                ole = getOutputs(inputValues)
                outputList.append(ole)
            except:
                dffix = int(math.floor((i-1)/2))
                if not centralDiff[dffix]:
                    # force it to the reported derivative.  This isn't great, but appears to be infrequent
                    # print('forcing derivative to mses value')
                    outputValues = {}
                    ks = ['camber','maxloc','thickness','Mach','Re',None]
                    if Alpha is not None:
                        ks[-1] = 'Alpha'
                        outputValues["CL"]            = outputList[0][0]['CL'] + outputList[0][1]['CL'][ks[dffix]] * deltaDict[ks[dffix]]
                    if CL is not None:
                        ks[-1] = 'CL'
                        outputValues["Alpha"]         = outputList[0][0]['Alpha'] + outputList[0][1]['Alpha'][ks[dffix]] * deltaDict[ks[dffix]]
                    outputValues["CD"]                = outputList[0][0]['CD']    + outputList[0][1]['CD'][ks[dffix]]    * deltaDict[ks[dffix]]
                    outputValues["CD_p"]              = outputList[0][0]['CD_p']  + outputList[0][1]['CD_p'][ks[dffix]]  * deltaDict[ks[dffix]]
                    outputValues["CD_v"]              = outputList[0][0]['CD_v']  + outputList[0][1]['CD_v'][ks[dffix]]  * deltaDict[ks[dffix]]
                    outputValues["CD_w"]              = outputList[0][0]['CD_w']  + outputList[0][1]['CD_w'][ks[dffix]]  * deltaDict[ks[dffix]]
                    outputValues["CM"]                = outputList[0][0]['CM']    + outputList[0][1]['CM'][ks[dffix]]    * deltaDict[ks[dffix]]
                    outputList.append([outputValues,{}])
                else:
                    if len(outputList) > 0:
                        outputList.append(outputList[0])
                    else:
                        raise RuntimeError('MSES did not converge')
                centralDiff[dffix] = False    
        else:
            outputList.append(outputList[0])


    delta_camber    += centralDiff[0] * delta_camber
    delta_maxloc    += centralDiff[1] * delta_maxloc    
    delta_thickness += centralDiff[2] * delta_thickness 
    delta_mach      += centralDiff[3] * delta_mach      
    delta_re        += centralDiff[4] * delta_re  
    if Alpha is not None:
        delta_alpha += centralDiff[5] * delta_alpha     
    if CL is not None:
        delta_cl    += centralDiff[5] * delta_cl

    outputValues = {}
    outputDerivatives = {}

    if Alpha is not None:
        outputValues["CL"]            = outputList[0][0]["CL"]
    if CL is not None:
        outputValues["Alpha"]         = outputList[0][0]["Alpha"]
        
    outputValues["CD"]                = outputList[0][0]["CD"]
    outputValues["CD_p"]              = outputList[0][0]["CD_p"]
    outputValues["CD_v"]              = outputList[0][0]["CD_v"]
    outputValues["CD_w"]              = outputList[0][0]["CD_w"]
    outputValues["CM"]                = outputList[0][0]["CM"]

    camber_ixl    = 1
    camber_ixu    = 2
    maxloc_ixl    = 3
    maxloc_ixu    = 4
    thickness_ixl = 5
    thickness_ixu = 6
    Mach_ixl      = 7
    Mach_ixu      = 8
    Re_ixl        = 9
    Re_ixu        = 10
    AC_ixl        = 11
    AC_ixu        = 12

    outputDerivatives["CD"] = {}
    outputDerivatives["CD"]["Mach"]      = (outputList[Mach_ixu][0]['CD']      - outputList[Mach_ixl][0]['CD']     ) / delta_mach
    outputDerivatives["CD"]["Re"]        = (outputList[Re_ixu][0]['CD']        - outputList[Re_ixl][0]['CD']       ) / delta_re
    outputDerivatives["CD"]["camber"]    = (outputList[camber_ixu][0]['CD']    - outputList[camber_ixl][0]['CD']   ) / delta_camber
    outputDerivatives["CD"]["maxloc"]    = (outputList[maxloc_ixu][0]['CD']    - outputList[maxloc_ixl][0]['CD']   ) / delta_maxloc
    outputDerivatives["CD"]["thickness"] = (outputList[thickness_ixu][0]['CD'] - outputList[thickness_ixl][0]['CD']) / delta_thickness

    outputDerivatives["CD_p"] = {}
    outputDerivatives["CD_p"]["Mach"]      = (outputList[Mach_ixu][0]['CD_p']      - outputList[Mach_ixl][0]['CD_p']     ) / delta_mach
    outputDerivatives["CD_p"]["Re"]        = (outputList[Re_ixu][0]['CD_p']        - outputList[Re_ixl][0]['CD_p']       ) / delta_re
    outputDerivatives["CD_p"]["camber"]    = (outputList[camber_ixu][0]['CD_p']    - outputList[camber_ixl][0]['CD_p']   ) / delta_camber
    outputDerivatives["CD_p"]["maxloc"]    = (outputList[maxloc_ixu][0]['CD_p']    - outputList[maxloc_ixl][0]['CD_p']   ) / delta_maxloc
    outputDerivatives["CD_p"]["thickness"] = (outputList[thickness_ixu][0]['CD_p'] - outputList[thickness_ixl][0]['CD_p']) / delta_thickness

    outputDerivatives["CD_v"] = {}
    outputDerivatives["CD_v"]["Mach"]      = (outputList[Mach_ixu][0]['CD_v']      - outputList[Mach_ixl][0]['CD_v']     ) / delta_mach
    outputDerivatives["CD_v"]["Re"]        = (outputList[Re_ixu][0]['CD_v']        - outputList[Re_ixl][0]['CD_v']       ) / delta_re
    outputDerivatives["CD_v"]["camber"]    = (outputList[camber_ixu][0]['CD_v']    - outputList[camber_ixl][0]['CD_v']   ) / delta_camber
    outputDerivatives["CD_v"]["maxloc"]    = (outputList[maxloc_ixu][0]['CD_v']    - outputList[maxloc_ixl][0]['CD_v']   ) / delta_maxloc
    outputDerivatives["CD_v"]["thickness"] = (outputList[thickness_ixu][0]['CD_v'] - outputList[thickness_ixl][0]['CD_v']) / delta_thickness

    outputDerivatives["CD_w"] = {}
    outputDerivatives["CD_w"]["Mach"]      = (outputList[Mach_ixu][0]['CD_w']      - outputList[Mach_ixl][0]['CD_w']     ) / delta_mach
    outputDerivatives["CD_w"]["Re"]        = (outputList[Re_ixu][0]['CD_w']        - outputList[Re_ixl][0]['CD_w']       ) / delta_re
    outputDerivatives["CD_w"]["camber"]    = (outputList[camber_ixu][0]['CD_w']    - outputList[camber_ixl][0]['CD_w']   ) / delta_camber
    outputDerivatives["CD_w"]["maxloc"]    = (outputList[maxloc_ixu][0]['CD_w']    - outputList[maxloc_ixl][0]['CD_w']   ) / delta_maxloc
    outputDerivatives["CD_w"]["thickness"] = (outputList[thickness_ixu][0]['CD_w'] - outputList[thickness_ixl][0]['CD_w']) / delta_thickness

    outputDerivatives["CM"] = {}
    outputDerivatives["CM"]["Mach"]      = (outputList[Mach_ixu][0]['CM']      - outputList[Mach_ixl][0]['CM']     ) / delta_mach
    outputDerivatives["CM"]["Re"]        = (outputList[Re_ixu][0]['CM']        - outputList[Re_ixl][0]['CM']       ) / delta_re
    outputDerivatives["CM"]["camber"]    = (outputList[camber_ixu][0]['CM']    - outputList[camber_ixl][0]['CM']   ) / delta_camber
    outputDerivatives["CM"]["maxloc"]    = (outputList[maxloc_ixu][0]['CM']    - outputList[maxloc_ixl][0]['CM']   ) / delta_maxloc
    outputDerivatives["CM"]["thickness"] = (outputList[thickness_ixu][0]['CM'] - outputList[thickness_ixl][0]['CM']) / delta_thickness

    if Alpha is not None:
        outputDerivatives["CL"] = {}
        outputDerivatives["CL"]["Alpha"]     = (outputList[AC_ixu][0]['CL']        - outputList[AC_ixl][0]['CL']       ) / delta_alpha
        outputDerivatives["CL"]["Mach"]      = (outputList[Mach_ixu][0]['CL']      - outputList[Mach_ixl][0]['CL']     ) / delta_mach
        outputDerivatives["CL"]["Re"]        = (outputList[Re_ixu][0]['CL']        - outputList[Re_ixl][0]['CL']       ) / delta_re
        outputDerivatives["CL"]["camber"]    = (outputList[camber_ixu][0]['CL']    - outputList[camber_ixl][0]['CL']   ) / delta_camber
        outputDerivatives["CL"]["maxloc"]    = (outputList[maxloc_ixu][0]['CL']    - outputList[maxloc_ixl][0]['CL']   ) / delta_maxloc
        outputDerivatives["CL"]["thickness"] = (outputList[thickness_ixu][0]['CL'] - outputList[thickness_ixl][0]['CL']) / delta_thickness
        
        outputDerivatives["CD"]["Alpha"]     = (outputList[AC_ixu][0]['CD']        - outputList[AC_ixl][0]['CD']       ) / delta_alpha
        outputDerivatives["CD_p"]["Alpha"]   = (outputList[AC_ixu][0]['CD_p']      - outputList[AC_ixl][0]['CD_p']     ) / delta_alpha
        outputDerivatives["CD_v"]["Alpha"]   = (outputList[AC_ixu][0]['CD_v']      - outputList[AC_ixl][0]['CD_v']     ) / delta_alpha
        outputDerivatives["CD_w"]["Alpha"]   = (outputList[AC_ixu][0]['CD_w']      - outputList[AC_ixl][0]['CD_w']     ) / delta_alpha
        outputDerivatives["CM"]["Alpha"]     = (outputList[AC_ixu][0]['CM']        - outputList[AC_ixl][0]['CM']       ) / delta_alpha
        
    if CL is not None:
        outputDerivatives["Alpha"] = {}
        outputDerivatives["Alpha"]["CL"]        = (outputList[AC_ixu][0]['Alpha']        - outputList[AC_ixl][0]['Alpha']       ) / delta_cl
        outputDerivatives["Alpha"]["Mach"]      = (outputList[Mach_ixu][0]['Alpha']      - outputList[Mach_ixl][0]['Alpha']     ) / delta_mach
        outputDerivatives["Alpha"]["Re"]        = (outputList[Re_ixu][0]['Alpha']        - outputList[Re_ixl][0]['Alpha']       ) / delta_re
        outputDerivatives["Alpha"]["camber"]    = (outputList[camber_ixu][0]['Alpha']    - outputList[camber_ixl][0]['Alpha']   ) / delta_camber
        outputDerivatives["Alpha"]["maxloc"]    = (outputList[maxloc_ixu][0]['Alpha']    - outputList[maxloc_ixl][0]['Alpha']   ) / delta_maxloc
        outputDerivatives["Alpha"]["thickness"] = (outputList[thickness_ixu][0]['Alpha'] - outputList[thickness_ixl][0]['Alpha']) / delta_thickness
        
        outputDerivatives["CD"]["CL"]     = (outputList[AC_ixu][0]['CD']        - outputList[AC_ixl][0]['CD']       ) / delta_cl
        outputDerivatives["CD_p"]["CL"]   = (outputList[AC_ixu][0]['CD_p']      - outputList[AC_ixl][0]['CD_p']     ) / delta_cl
        outputDerivatives["CD_v"]["CL"]   = (outputList[AC_ixu][0]['CD_v']      - outputList[AC_ixl][0]['CD_v']     ) / delta_cl
        outputDerivatives["CD_w"]["CL"]   = (outputList[AC_ixu][0]['CD_w']      - outputList[AC_ixl][0]['CD_w']     ) / delta_cl
        outputDerivatives["CM"]["CL"]     = (outputList[AC_ixu][0]['CM']        - outputList[AC_ixl][0]['CM']       ) / delta_cl
    
    myProblem.close()
    
    # sensxFile = os.path.join(workDir,'Scratch','msesAIM0','sensx.airfoil')
    # sensDict = readSensx(sensxFile)

    if cleanup:
        if os.path.exists(workDir) and os.path.isdir(workDir):
            shutil.rmtree(workDir)
        
    return [outputValues, outputDerivatives, outputList[0][1]]#, sensDict]

# import matplotlib.pyplot as plt
# import numpy as np
      
# sweepParamList = ["Mach","Re","Alpha","camber","maxloc","thickness"]

# for sweepParam in sweepParamList:
#     Mach = 0.3
#     Re = 1e7
#     Alpha = 1.0
#     camber = 0.02
#     maxloc = 0.40
#     thickness = 0.12

#     pltTitle = "Plots varying %s from a baseline NACA %d%d%d at %.1f deg AOA, Mach %.2f, and Re %.1e"%(sweepParam.lower(),int(camber*100),int(maxloc*10),int(thickness*100),Alpha,Mach,Re)

#     if sweepParam == "Mach":
#         vals = np.linspace(0.2,0.6,9)
#     elif sweepParam == "Re":
#         vals = 10**np.linspace(6,7.5,11)
#     elif sweepParam == "Alpha":
#         vals = np.linspace(0,10,13)
#     elif sweepParam == "camber":
#         vals = np.linspace(0,.06,7)
#     elif sweepParam == "maxloc":
#         vals = np.linspace(.1,.6,6)
#     elif sweepParam == "thickness":
#         vals = np.linspace(0.07,0.25,10)
#     else:
#         cody
        
#     outputData = []
#     for vl in vals:
#         if sweepParam == "Mach":
#             Mach = vl
#         elif sweepParam == "Re":
#             Re = vl
#         elif sweepParam == "Alpha":
#             Alpha = vl
#         elif sweepParam == "camber":
#             camber = vl
#         elif sweepParam == "maxloc":
#             maxloc = vl
#         elif sweepParam == "thickness":
#             thickness = vl
#         else:
#             cody

#         opt = msesRunNACA4_Alpha(camber, maxloc, thickness, Mach, Re, Alpha)
#         outputData.append(opt)

#     fig , ax = plt.subplots(figsize=(24,12), constrained_layout=True)
#     plt.rcParams['font.size'] = 18
#     fig.suptitle(pltTitle)
#     ylabs = ['CL','CD','CM']
#     for i in range(0,len(ylabs)):
#         plt.subplot(2,3,i+1)
#         resVals = np.array([rw[0][ylabs[i]] for rw in outputData])
#         plt.plot(vals,resVals)
#         plt.xlabel("%s"%(sweepParam.lower()))
#         plt.ylabel(ylabs[i])
#         plt.grid(1)

#         plt.subplot(2,3,i+4)
#         derApp = (resVals[1:]-resVals[0:-1])/(vals[1]-vals[0])
#         vals_avg = (vals[1:]+vals[0:-1])/2
#         plt.plot(vals_avg,derApp)
#         plt.plot(vals, np.array([rw[1][ylabs[i]][sweepParam] for rw in outputData]))

#         plt.xlabel("%s"%(sweepParam.lower()))
#         plt.ylabel("d(%s)/d(%s)"%(ylabs[i],sweepParam.lower()))
#         plt.legend(['FD Approx','CAPS Reported'])
#         plt.grid(1)

#     plt.savefig("NACA4_%s_Study.png"%(sweepParam))
