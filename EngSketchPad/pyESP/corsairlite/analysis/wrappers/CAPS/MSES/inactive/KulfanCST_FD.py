import pyCAPS
import os
import numpy as np
import shutil
import math
import copy
from corsairlite.analysis.wrappers.CAPS.MSES.readSensx import readSensx

def msesRunKulfan_Alpha(Aupper,Alower,Mach,Re,Alpha,
                         Acrit = None, GridAlpha = None, 
                         Coarse_Iteration = None, Fine_Iteration = None,
                         xGridRange = None, yGridRange = None, cleanup = True ):
    
    return msesRunKulfan( Aupper,Alower,
                         Mach, Re,        
                         Alpha = Alpha,
                         CL = None,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )
    
def msesRunKulfan_CL(Aupper,Alower,Mach,Re,CL,
                         Acrit = None, GridAlpha = None, 
                         Coarse_Iteration = None, Fine_Iteration = None,
                         xGridRange = None, yGridRange = None, cleanup = True ):
    
    return msesRunKulfan( Aupper,Alower,
                         Mach, Re,        
                         Alpha = None,
                         CL = CL,
                         Acrit = Acrit, GridAlpha = GridAlpha, 
                         Coarse_Iteration = Coarse_Iteration, Fine_Iteration = Fine_Iteration,
                         xGridRange = xGridRange, yGridRange = yGridRange, cleanup = cleanup )


def msesRunKulfan(Aupper,Alower,
                 Mach, Re,        
                 Alpha = None,
                 CL = None,
                 Acrit = None, GridAlpha = None, 
                 Coarse_Iteration = None, Fine_Iteration = None,
                 xGridRange = None, yGridRange = None, cleanup = True ):

    geometryScript = "Kulfan_airfoilSection.csm"
    runDirectory = "msesAnalysisTest"

    pstr = ""
    pstr += 'attribute capsAIM $xfoilAIM;tsfoilAIM;msesAIM \n'
    pstr += 'dimension class 1 2 \n'
    pstr += 'dimension ztail 1 2 \n'
    pstr += 'dimension aupper 1 %d \n'%(len(Aupper))
    pstr += 'dimension alower 1 %d \n'%(len(Alower))
    pstr += 'despmtr class "0.5; 1.0;" \n'
    pstr += 'despmtr ztail "0.001; -0.001;" \n'
    pstr += 'despmtr aupper " '
    for vl in Aupper:
        pstr += " %f;"%(vl)
    pstr += '" \n'
    pstr += 'despmtr alower " '
    for vl in Alower:
        pstr += " %f;"%(vl)
    pstr += '" \n'
    pstr += 'udparg kulfan class class \n'
    pstr += 'udparg kulfan ztail ztail \n'
    pstr += 'udparg kulfan aupper aupper \n'
    pstr += 'udparg kulfan alower alower \n'
    pstr += 'udprim kulfan \n'
    pstr += 'extract 0\n'

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
        mses.input.GridAlpha        = float(math.ceil(180.0 * CL / (4*np.pi**2))) # Approximating alpha wiht thin airfoil theory, best thought I have
    if Coarse_Iteration is not None:
        mses.input.Coarse_Iteration = Coarse_Iteration # Default 20
    if Fine_Iteration is not None:
        mses.input.Fine_Iteration   = Fine_Iteration   # Default 200
    if xGridRange is not None:
        mses.input.xGridRange       = xGridRange       # Default [-1.75, 2.75]
    if yGridRange is not None:
        mses.input.yGridRange       = yGridRange       # Default [-2.0, 2.5]

    # mses.input.Airfoil_Points = 231 # points

    if Mach <= 0.1 :
        mses.input.ISMOM = 2 # see mses documentation for M<=0.1

    mses.input.Design_Variable = {"aupper": {}, "alower":{}}

    delta_A         = 0.005
    delta_mach      = 0.0015
    delta_re        = 80
    delta_alpha     = 0.07
    delta_cl        = 0.01
    deltaDict = {}
    deltaDict['A']         = delta_A
    deltaDict['Mach']      = delta_mach
    deltaDict['Re']        = delta_re
    deltaDict['Alpha']     = delta_alpha
    deltaDict['CL']        = delta_cl
    centralDiff = [True]*len(Aupper) + [True]*len(Alower) + [True]*3


    inputList = []
    if Alpha is not None:
        inputStack0 = np.array(Aupper).tolist() + np.array(Alower).tolist() + [Mach, Re, Alpha]
        inputList.append(inputStack0)
    if CL is not None:
        inputStack0 = np.array(Aupper).tolist() + np.array(Alower).tolist() + [Mach, Re, CL]
        inputList.append(inputStack0)
    inputStack = copy.deepcopy(inputStack0)


    # Kulfan Coefficients
    Astack = np.array(Aupper).tolist() + np.array(Alower).tolist()
    isIx = -1
    for i in range(0,len(Astack)):
        isIx += 1
        A_l = Astack[i] - delta_A
        A_u = Astack[i] + delta_A
        inputStack = copy.deepcopy(inputStack0)
        inputStack[isIx] = A_l
        inputList.append(inputStack)
        inputStack = copy.deepcopy(inputStack0)
        inputStack[isIx] = A_u
        inputList.append(inputStack)

    # Mach
    isIx += 1
    mach_l = Mach - delta_mach
    mach_u = Mach + delta_mach
    if mach_l < 0.0:
        inputList.append(None)
        centralDiff[isIx] = False
    else:
        inputStack = copy.deepcopy(inputStack0)
        inputStack[isIx] = mach_l
        inputList.append(inputStack)
    inputStack = copy.deepcopy(inputStack0)
    inputStack[isIx] = mach_u
    inputList.append(inputStack)

    # Re
    isIx += 1
    re_l = Re - delta_re
    re_u = Re + delta_re
    inputStack = copy.deepcopy(inputStack0)
    inputStack[isIx] = re_l
    inputList.append(inputStack)
    inputStack = copy.deepcopy(inputStack0)
    inputStack[isIx] = re_u
    inputList.append(inputStack)

    isIx += 1
    if Alpha is not None:
        alpha_l = Alpha - delta_alpha
        alpha_u = Alpha + delta_alpha
        inputStack = copy.deepcopy(inputStack0)
        inputStack[isIx] = alpha_l
        inputList.append(inputStack)
        inputStack = copy.deepcopy(inputStack0)
        inputStack[isIx] = alpha_u
        inputList.append(inputStack)
    if CL is not None:
        cl_l = CL - delta_cl
        cl_u = CL + delta_cl
        inputStack = copy.deepcopy(inputStack0)
        inputStack[isIx] = cl_l
        inputList.append(inputStack)
        inputStack = copy.deepcopy(inputStack0)
        inputStack[isIx] = cl_u
        inputList.append(inputStack)

    outputList = []
    for i in range(0,len(inputList)):
        inputValues = inputList[i]
        # print(inputValues)
        if inputValues is not None:
            try:
                Aupper_input = inputValues[0:len(Aupper)]
                Alower_input = inputValues[len(Aupper):len(Aupper) + len(Alower)]
                settingIx = len(Aupper) + len(Alower)
                mses.geometry.despmtr.aupper            = Aupper_input
                mses.geometry.despmtr.alower            = Alower_input
                mses.input.Mach              = inputValues[settingIx]
                mses.input.Re                = inputValues[settingIx + 1]

                if Alpha is not None:
                    mses.input.Alpha = inputValues[settingIx + 2]
                if CL is not None:
                    mses.input.CL = inputValues[settingIx + 2]

                firstOutput = mses.output.CD # this runs mses
                optDict = {}
                dervDict = {}
                msesKeys = mses.output.keys()
                for ky in msesKeys:
                    optDict[ky]=mses.output[ky].value
                    dervDict[ky] = mses.output[ky].deriv()

                outputList.append([optDict,dervDict])
            except: #Mses failed to converge.  Need to try left/right difference, then revert
                dffix = int(math.floor((i-1)/2))
                if not centralDiff[dffix]: # There has already been a failure causing central diff to flag off, so need to revert
                    # force it to the reported derivative.  This isn't great, but appears to be infrequent
                    # print('forcing derivative to mses value')
                    outputValues = {}
                    ks = ['aupper']*len(Aupper) + ['alower']*len(Alower) + ['Mach','Re',None]
                    if dffix > (len(Aupper) + len(Alower)):
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
                        if dffix < len(Aupper):
                            aKey = 'aupper'
                            aIx = dffix
                        else:
                            aKey = 'alower'
                            aIx = dffix-len(Aupper)
                            
                        if Alpha is not None:
                            ks[-1] = 'Alpha'
                            outputValues["CL"]            = outputList[0][0]['CL'] + outputList[0][1]['CL'][aKey][aIx] * deltaDict['A']
                        if CL is not None:
                            ks[-1] = 'CL'
                            outputValues["Alpha"]         = outputList[0][0]['Alpha'] + outputList[0][1]['Alpha'][aKey][aIx] * deltaDict['A']
                        outputValues["CD"]                = outputList[0][0]['CD']    + outputList[0][1]['CD'][aKey][aIx]    * deltaDict['A']
                        outputValues["CD_p"]              = outputList[0][0]['CD_p']  + outputList[0][1]['CD_p'][aKey][aIx]  * deltaDict['A']
                        outputValues["CD_v"]              = outputList[0][0]['CD_v']  + outputList[0][1]['CD_v'][aKey][aIx]  * deltaDict['A']
                        outputValues["CD_w"]              = outputList[0][0]['CD_w']  + outputList[0][1]['CD_w'][aKey][aIx]  * deltaDict['A']
                        outputValues["CM"]                = outputList[0][0]['CM']    + outputList[0][1]['CM'][aKey][aIx]    * deltaDict['A']
                        outputList.append([outputValues,{}])
                else: # one of the cases passed, so just need to turn off central diff and proceed
                    if len(outputList) > 0:
                        outputList.append(outputList[0])
                    else: # this is the actual evaluation point, return error
                        raise RuntimeError('MSES did not converge')
                centralDiff[dffix] = False 
        else: # someting like M<0 happened, so just append and move on
            outputList.append(outputList[0])

    centralDiff_Aupper = centralDiff[0:len(Aupper)]
    centralDiff_Alower = centralDiff[len(Aupper):len(Aupper)+len(Alower)]
    centralDiff_other  = centralDiff[len(Aupper)+len(Alower):]

    delta_aupper = [delta_A + dfv * delta_A for dfv in centralDiff_Aupper]
    delta_alower = [delta_A + dfv * delta_A for dfv in centralDiff_Alower]
    delta_mach      += centralDiff_other[0] * delta_mach      
    delta_re        += centralDiff_other[1] * delta_re  
    if Alpha is not None:
        delta_alpha += centralDiff_other[2] * delta_alpha     
    if CL is not None:
        delta_cl    += centralDiff_other[2] * delta_cl


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

    Mach_ixl      = 1 + 2*len(Aupper) + 2*len(Alower)
    Mach_ixu      = Mach_ixl + 1
    Re_ixl        = Mach_ixu + 1
    Re_ixu        = Re_ixl + 1
    AC_ixl        = Re_ixu + 1
    AC_ixu        = AC_ixl + 1

    for ky in ['CD','CD_p','CD_v','CD_w','CM']:
        outputDerivatives[ky] = {}
        outputDerivatives[ky]["Mach"]      = (outputList[Mach_ixu][0][ky]      - outputList[Mach_ixl][0][ky]     ) / delta_mach
        outputDerivatives[ky]["Re"]        = (outputList[Re_ixu][0][ky]        - outputList[Re_ixl][0][ky]       ) / delta_re
        aupper_dout = [None]*len(Aupper)
        for i in range(0,len(Aupper)):
            aupper_ixl = 1 + 2*i
            aupper_ixu = 1 + 2*i + 1
            aupper_dout[i] = (outputList[aupper_ixu][0][ky]        - outputList[aupper_ixl][0][ky]       ) / delta_aupper[i]
        outputDerivatives[ky]["aupper"] = np.array(aupper_dout)

        alower_dout = [None]*len(Alower)
        for i in range(0,len(Alower)):
            alower_ixl = 1 + 2*i +     2*len(Aupper)
            alower_ixu = 1 + 2*i + 1 + 2*len(Aupper)
            alower_dout[i] = (outputList[alower_ixu][0][ky]        - outputList[alower_ixl][0][ky]       ) / delta_alower[i]
        outputDerivatives[ky]["alower"] = np.array(alower_dout)
        
    if Alpha is not None:
        outputDerivatives["CL"] = {}
        outputDerivatives["CL"]["Alpha"]     = (outputList[AC_ixu][0]['CL']        - outputList[AC_ixl][0]['CL']       ) / delta_alpha
        outputDerivatives["CL"]["Mach"]      = (outputList[Mach_ixu][0]['CL']      - outputList[Mach_ixl][0]['CL']     ) / delta_mach
        outputDerivatives["CL"]["Re"]        = (outputList[Re_ixu][0]['CL']        - outputList[Re_ixl][0]['CL']       ) / delta_re
        aupper_dout = [None]*len(Aupper)
        for i in range(0,len(Aupper)):
            aupper_ixl = 1 + 2*i
            aupper_ixu = 1 + 2*i + 1
            aupper_dout[i] = (outputList[aupper_ixu][0]["CL"]        - outputList[aupper_ixl][0]["CL"]       ) / delta_aupper[i]
        outputDerivatives["CL"]["aupper"] = np.array(aupper_dout)
        alower_dout = [None]*len(Alower)
        for i in range(0,len(Alower)):
            alower_ixl = 1 + 2*i +     2*len(Aupper)
            alower_ixu = 1 + 2*i + 1 + 2*len(Aupper)
            alower_dout[i] = (outputList[alower_ixu][0]["CL"]        - outputList[alower_ixl][0]["CL"]       ) / delta_alower[i]
        outputDerivatives["CL"]["alower"] = np.array(alower_dout)

        for ky in ['CD','CD_p','CD_v','CD_w','CM']:
            outputDerivatives[ky]["Alpha"]     = (outputList[AC_ixu][0][ky]        - outputList[AC_ixl][0][ky]       ) / delta_alpha

    if CL is not None:
        outputDerivatives["Alpha"] = {}
        outputDerivatives["Alpha"]["CL"]        = (outputList[AC_ixu][0]['Alpha']        - outputList[AC_ixl][0]['Alpha']       ) / delta_cl
        outputDerivatives["Alpha"]["Mach"]      = (outputList[Mach_ixu][0]['Alpha']      - outputList[Mach_ixl][0]['Alpha']     ) / delta_mach
        outputDerivatives["Alpha"]["Re"]        = (outputList[Re_ixu][0]['Alpha']        - outputList[Re_ixl][0]['Alpha']       ) / delta_re
        aupper_dout = [None]*len(Aupper)
        for i in range(0,len(Aupper)):
            aupper_ixl = 1 + 2*i
            aupper_ixu = 1 + 2*i + 1
            aupper_dout[i] = (outputList[aupper_ixu][0]["Alpha"]        - outputList[aupper_ixl][0]["Alpha"]       ) / delta_aupper[i]
        outputDerivatives["Alpha"]["aupper"] = np.array(aupper_dout)
        alower_dout = [None]*len(Alower)
        for i in range(0,len(Alower)):
            alower_ixl = 1 + 2*i +     2*len(Aupper)
            alower_ixu = 1 + 2*i + 1 + 2*len(Aupper)
            alower_dout[i] = (outputList[alower_ixu][0]["Alpha"]        - outputList[alower_ixl][0]["Alpha"]       ) / delta_alower[i]
        outputDerivatives["Alpha"]["alower"] = np.array(alower_dout)

        for ky in ['CD','CD_p','CD_v','CD_w','CM']:
            outputDerivatives[ky]["CL"]     = (outputList[AC_ixu][0][ky]        - outputList[AC_ixl][0][ky]       ) / delta_cl

    myProblem.close()
    
    # sensxFile = os.path.join(workDir,'Scratch','msesAIM0','sensx.airfoil')
    # sensDict = readSensx(sensxFile)

    if cleanup:
        if os.path.exists(workDir) and os.path.isdir(workDir):
            shutil.rmtree(workDir)

    return [outputValues, outputDerivatives]#, sensDict]

