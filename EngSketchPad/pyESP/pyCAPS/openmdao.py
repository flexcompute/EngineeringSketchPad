from . import caps
import os

def getPythonObjShape(dataValue):
    
    numRow = 0
    numCol = 0
    
    if not isinstance(dataValue, list):
        dataValue = [dataValue] # Convert to list
            
    numRow = len(dataValue)
        
    if isinstance(dataValue[0], list):
        numCol = len(dataValue[0]) 
    else:
        numCol = 1
        
    # Make sure shape is consistent
    if numCol != 1:
        for i in range(numRow):
            if len(dataValue[i]) != numCol:
                print("Inconsistent list sizes!", i, dataValue[i], numCol, dataValue[0])
                raise CAPSError(cCAPS.CAPS_BADVALUE, msg="Inconsistent list sizes! {!r} {!r} {!r} {!r}".format(i, dataValue[i], numCol, dataValue[0]))
    
    return numRow, numCol


def setElementofCAPSTuple(tupleValue, tupleKey, elementValue, dictKey = None):
    
    if not isinstance(tupleValue, list):
        tupleValue = [tupleValue]
        
    output = []
    for tupleElement in tupleValue:
        
        key = tupleElement[0]
        value = tupleElement[1]
        
        if tupleKey == key:
            
            if isinstance(value, dict): # Modify element of a dictionary 
                value[dictKey] = elementValue
            else: # Else just modify the value 
                value = elementValue        
        
        output.append((key, value))
    
    return output


def createOpenMDAOComponent(analysis, 
                            inputParam, 
                            outputParam, 
                            executeCommand = None, 
                            inputFile = None, 
                            outputFile = None, 
                            sensitivity = {}, 
                            changeDir = None, 
                            saveIteration = False):
    
    try: 
        from openmdao import __version__ as version
    except:
        raise caps.CAPSError(caps.CAPS_NOTFOUND, "Error: Unable to determine OpenMDAO version!")
    try:
        from numpy import ndarray
    except:
        raise caps.CAPSError(caps.CAPS_NOTFOUND, "Error: Unable to import numpy!")
    
    print("OpenMDAO version - ", version)
    versionMajor = int(version.split(".")[0])
    versionMinor = int(version.split(".")[1])
    
    if versionMajor >=2:
        try: 
            from  openmdao.api import ExternalCodeComp, ExplicitComponent     
        except: 
            raise caps.CAPSError(caps.CAPS_NOTFOUND, "Error: Unable to import ExternalCodeComp or ExplicitComponent from openmdao.api!")
    
    elif versionMajor == 1 and versionMinor == 7:   
        try: 
            from  openmdao.api import ExternalCode, Component     
        except: 
            raise caps.CAPSError(caps.CAPS_NOTFOUND, "Error: Unable to import ExternalCode or Component from openmdao.api!")
    else:
        raise caps.CAPSError(caps.CAPS_BADVALUE, "Error: Unsupported OpenMDAO version!")
    
    # Change to lists if not lists already 
    if not isinstance(inputParam, list):
        inputParam = [inputParam]
        
    if not isinstance(outputParam, list):
        outputParam = [outputParam]

    def filterInputs(varName, validInput, validGeometry):
        
        value = None
        
        if ":" in varName: # Is the input trying to modify a Tuple?
          
            analysisTuple = varName.split(":")[0] 
            if analysisTuple not in validInput.keys():
                if varName not in validGeometry.keys(): #Check geometry
                    raise caps.CAPSError(caps.CAPS_BADVALUE, "It appears the input parameter", varName, "is trying to modify a Tuple input, but", analysisTuple,
                          "this is not a valid ANALYSISIN or GEOMETRYIN variable! It will not be added to the OpenMDAO component")
                else: 
                    value = validGeometry[varName]
            else:  
                
                if ":" in varName:
                    tuples = validInput[analysisTuple].value
                    if not isinstance(tuples,list):
                        tuples = [tuples]
                    
                    tupleValue = varName.split(":")[1] 
                    
                    for i in tuples:
                        if tupleValue == i[0]:
                            if varName.count(":") == 2:
                                dictValue = varName.split(":")[2] 
                                value = i[1][dictValue]
                            elif varName.count(":") == 1:
                                value = i[1]
        else:
            # Check to make sure variable is in analysis
            if varName not in validInput.keys():
                if varName not in validGeometry.keys():
                    raise caps.CAPSError(caps.CAPS_BADVALUE, "Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!", 
                          "It will not be added to the OpenMDAO component.")
                else:
                    value = validGeometry[varName].value
            else:
                value = validInput[varName].value
                
        # Check for None
        if value is None:
            value = 0.0
        
        return value
    
    def filterOutputs(varName, validOutput):
        
        value = None
        
        # Check to make sure variable is in analysis
        if varName not in validOutput:
            raise caps.CAPSError(caps.CAPS_BADVALUE, "Output parameter '" + str(varName) + "' is not a valid ANALYSISOUT variable!")
        
        # Check for None
        if value is None:
            value = 0.0
        
        return value
    
    def setInputs(varName, varValue, analysis, validInput, validGeometry):
       
        if isinstance(varValue, ndarray):
            varValue = varValue.tolist()
        
        if isinstance(varValue,list):
            numRow, numCol = getPythonObjShape(varValue)
            if numRow == 1 and numCol == 1:
                varValue = varValue[0]
            elif numRow == 1 and numCol > 1:
                varValue = varValue[0]
                
        if ":" in varName: # Is the input trying to modify a Analysis Tuple?
            
            analysisTuple = varName.split(":")[0] 
            
            if analysisTuple not in validInput:
                if varName in validGeometry: #Check geometry
                    analysis.geometry.despmtr[varName].value = varValue
                else:
                    print("Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!", 
                          "It will not be modified to the OpenMDAO component. This should have already been reported!")
            else:
                if varName.count(":") == 2:
                   
                    value = setElementofCAPSTuple(analysis.getAnalysisVal(analysisTuple), # Validity has already been checked
                                                  varName.split(":")[1], # Tuple key
                                                  varValue, # Value to set 
                                                  varName.split(":")[2]) # Dictionary key
                else: 
                    
                    value = setElementofCAPSTuple(analysis.getAnalysisVal(analysisTuple), # Validity has already been checked
                                                  varName.split(":")[1],  # Tuple key
                                                  varValue) # Value to set 
                
                analysis.input[analysisTuple].value = value
            
        else:
            
            if varName in validInput:
                
                analysis.input[varName].value = varValue
            
            elif varName in validGeometry:
                
                analysis.geometry.despmtr[varName].value = varValue
            else:
                raise caps.CAPSError(caps.CAPS_BADVALUE, "Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!", 
                                                          "It will not be modified to the OpenMDAO component. This should have already been reported!")    
    
#     def checkParentExecution(analysis, checkInfo = False):
#         
#         def infoDictCheck(analysis):
#             
#             infoDict = analysis.getAnalysisInfo(printInfo = False, infoDict = True)
#         
#             if infoDict["executionFlag"] == True:
#                 return True
#             else:
#                 print(analysis.aimName, "can't self execute!")
#                 return False
#         
#         if analysis.parents: # If we have parents 
#             
#             # Check analysis parents
#             for parent in analysis.parents:
#         
#                 if not checkParentExecution(analysis.capsProblem.analysis[parent], True):
#                     return False
#                 
#                 if checkInfo:
#                     return infoDictCheck(analysis)
#             
#         if checkInfo:
#             return infoDictCheck(analysis)
#     
    def executeLinked(analysis):
         
        links = analysis.dirty()
             
        for link in links:
            link.preAnalysis()
            link.postAnalysis()

        return 
    
    def initialInstance(analysisDir):
        
         # Get current directory
        currentDir = os.getcwd()
        
        # Change to analysis directory
        os.chdir(analysisDir)
        
        # Determine instance number
        instance = 0
        key = "Instance_"
        
        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]
        
        for i in currentFiles:
            
            path = os.path.join(os.getcwd(), i)

            if os.path.isdir(path) and key in i:
                keepFiles.remove(i)
                
                if int(i.replace(key, "")) >= instance:
                    instance = int(i.replace(key, "")) + 1

        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return 

        # Make new instance directory
        instance = os.path.join(os.getcwd(), key + str(instance))
        os.mkdir(instance)
        
        for i in keepFiles:
            # print"Change file - ", os.path.join(instance,i)
            os.rename(i, os.path.join(instance,i))
                
        # Change back to current directory
        os.chdir(currentDir)
        
    def storeIteration(analysisDir):
        
        # Get current directory
        currentDir = os.getcwd()
        
        # Change to analysis directory
        os.chdir(analysisDir)
        
        # Determine instance number
        interation = 0
        key = "Iteration_"
        key2 = "Instance_"
        
        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]
        
        for i in currentFiles:
            path = os.path.join(os.getcwd(), i)

            if os.path.isdir(path):
                
                if key in i or key2 in i:
                    keepFiles.remove(i)
                
                if key in i:
                    if int(i.replace(key, "")) >= interation:
                        interation = int(i.replace(key, "")) + 1

        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return 
        
        # Make new interation directory
        interation = os.path.join(os.getcwd(), key + str(interation))
        if not os.path.isdir(interation):
            os.mkdir(interation)
        
        # Move files into new interation folder
        for i in keepFiles:
            # print"Change file - ", os.path.join(interation,i)
            os.rename(i, os.path.join(interation,i))
        
        # Change back to current directory
        os.chdir(currentDir)
    
    if versionMajor >=2:
        class openMDAOExternal_V2(ExternalCodeComp):      
            def __init__(self, 
                         analysis, 
                         inputParam, outputParam, executeCommand, 
                         inputFile = None, outputFile = None, 
                         sensitivity = {}, 
                         changeDir = True,
                         saveIteration = False):
                
                super(openMDAOExternal_V2, self).__init__()
                
                self.analysis = analysis
                self.inputParam = inputParam
                self.outputParam = outputParam
                self.executeCommand = executeCommand
                self.inputFile = inputFile
                self.outputFile = outputFile
                self.sensitivity = sensitivity
                self.changeDir = changeDir
                self.saveIteration = saveIteration
                
                # Get check to make sure         
                self.validInput    = self.analysis.input
                self.validOutput   = self.analysis.output
                self.validGeometry = self.analysis.geometry.despmtr
            
            def setup(self):
                
                for i in self.inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
             
                    self.add_input(i, val=value)
                
                for i in self.outputParam:
                    value = filterOutputs(i, self.validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                    
                if self.inputFile is not None:
                    if isinstance(self.inputFile, list):
                        self.options['external_input_files'] = self.inputFile
                    else:
                        self.options['external_input_files'] = [self.inputFile]
                
                if self.outputFile is not None:
                    if isinstance(self.outputFile, list):
                        self.options['external_output_files'] = self.outputFile
                    else:
                        self.options['external_output_files'] = [self.outputFile]
               
                if isinstance(self.executeCommand, list):
                    self.options['command'] = self.executeCommand
                else:
                    self.options['command'] = [self.executeCommand]
                    
                method = "fd"
                step= None
                form = None
                step_calc = None
                
                for i in self.sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = self.sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = self.sensitivity[i]
                    elif "form" in i: 
                        form = self.sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = self.sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)
                  
                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return cCAPS.CAPS_BADVALUE
                
#                  self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                self.declare_partials('*', '*', method=method, step=step, form=form, step_calc=step_calc)
                
                #self.changeDir = changeDir
                
                #self.saveIteration = saveIteration     
                
                if self.saveIteration:
                    initialInstance(self.analysis.analysisDir)           
                     
            def compute(self, params, unknowns):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysis, self.validInput, self.validGeometry)
                
                # Run parents 
                executeLinked(self.analysis)
                
                if self.saveIteration:
                    storeIteration(self.analysis.analysisDir)
                        
                self.analysis.preAnalysis()
                    
                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise 
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysis.analysisDir)
                    
                #Parent compute function actually runs the external code
                super(openMDAOExternal_V2, self).compute(params, unknowns)
                
                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None: 
                    os.chdir(currentDir)
                    
                self.analysis.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysis.getAnalysisOutVal(varname = i)
            
        class openMDAOExternal_V2ORIGINAL(ExternalCodeComp):      
            def __init__(self, 
                         analysis, 
                         inputParam, outputParam, executeCommand, 
                         inputFile = None, outputFile = None, 
                         sensitivity = {}, 
                         changeDir = True,
                         saveIteration = False):
                
                super(openMDAOExternal_V2, self).__init__()
                
                self.analysis = analysis
                
                # Get check to make sure         
                self.validInput  = self.analysis.analysis.input
                validOutput      = self.analysis.output
                self.validGeometry = self.analysis.geometry.despmtr
                
                for i in inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
             
                    self.add_input(i, val=value)
                
                for i in outputParam:
                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                    
                if inputFile is not None:
                    if isinstance(inputFile, list):
                        self.options['external_input_files'] = inputFile
                    else:
                        self.options['external_input_files'] = [inputFile]
                
                if outputFile is not None:
                    if isinstance(outputFile, list):
                        self.options['external_output_files'] = outputFile
                    else:
                        self.options['external_output_files'] = [outputFile]
               
                if isinstance(executeCommand, list):
                    self.options['command'] = executeCommand
                else:
                    self.options['command'] = [executeCommand]
                    
                method = "fd"
                step= None
                form = None
                step_calc = None
                
                for i in sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = sensitivity[i]
                    elif "form" in i: 
                        form = sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)
                  
                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return cCAPS.CAPS_BADVALUE
                
#                 self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                self.declare_partials('*', '*', method=method, step=step, form=form, step_calc=step_calc)
                
                self.changeDir = changeDir
                
                self.saveIteration = saveIteration     
                
                if self.saveIteration:
                    initialInstance(self.analysis.analysisDir)           
                     
            def compute(self, params, unknowns):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysis, self.validInput, self.validGeometry)
                
                # Run parents 
                executeLinked(self.analysis)
                
                if self.saveIteration:
                    storeIteration(self.analysis.analysisDir)
                        
                self.analysis.preAnalysis()
                    
                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise 
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysis.analysisDir)
                    
                #Parent compute function actually runs the external code
                super(openMDAOExternal_V2, self).compute(params, unknowns)
                
                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None: 
                    os.chdir(currentDir)
                    
                self.analysis.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysis.getAnalysisOutVal(varname = i)
            
            # TO add sensitivities       
    #         def compute_partials(self, inputs, partials):
    #             outputs = {}
    #      the parent compute function actually runs the external code
    #         super(ParaboloidExternalCodeCompDerivs, self).compute(inputs, outputs)
    # 
    #         # parse the derivs file from the external code and set partials
    #         with open(self.derivs_file, 'r') as derivs_file:
    #             partials['f_xy', 'x'] = float(derivs_file.readline())
    #             partials['f_xy', 'y'] = float(derivs_file.readline())
                
        class openMDAOComponent_V2(ExplicitComponent):
                
            def __init__(self, 
                         analysis, 
                         inputParam, outputParam, 
                         sensitivity = {}, 
                         saveIteration = False):
                    
                super(openMDAOComponent_V2, self).__init__()
    
                self.analysis = analysis
                    
                # Get check to make sure         
                self.validInput  = self.analysis.input
                validOutput      = self.analysis.output
                self.validGeometry = self.analysis.geometry.despmtr
                
                for i in inputParam:
                    
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
                    
                    self.add_input(i, val=value)
                
                for i in outputParam:
                    
                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                   
                method = "fd"
                step= None
                form = None
                step_calc = None
                
                for i in sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = sensitivity[i]
                    elif "form" in i: 
                        form = sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)
                
                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return cCAPS.CAPS_BADVALUE
                
                self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                
                self.saveIteration = saveIteration
                
                if self.saveIteration:
                    initialInstance(self.analysis.analysisDir)  
               
            def compute(self, params, unknowns):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysis, self.validInput, self.validGeometry)
               
                # Run parents 
                executeLinked(self.analysis)
                
                if self.saveIteration:
                    storeIteration(self.analysis.analysisDir)
                    
                self.analysis.preAnalysis()
                self.analysis.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysis.getAnalysisOutVal(varname = i)
            
            # TO add sensitivities       
    #         def compute_partials(self, inputs, partials):
    #             partials['area', 'length'] = inputs['width']
    #             partials['area', 'width'] = inputs['length']
    
    else:
        class openMDAOExternal(ExternalCode):      
            def __init__(self, 
                         analysis, 
                         inputParam, outputParam, executeCommand, 
                         inputFile = None, outputFile = None, 
                         sensitivity = {}, 
                         changeDir = True,
                         saveIteration = False):
                
                super(openMDAOExternal, self).__init__()
                
                self.analysis = analysis
                
                # Get check to make sure         
                self.validInput  = self.analysis.input
                validOutput      = self.analysis.output
                self.validGeometry = self.analysis.geometry.despmtr
                
                for i in inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
                    
                    self.add_param(i, val=value)
                
                for i in outputParam:
                    value = filterOutputs(i, validOutput)
                    
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                    
                if inputFile is not None:
                    if isinstance(inputFile, list):
                        self.options['external_input_files'] = inputFile
                    else:
                        self.options['external_input_files'] = [inputFile]
                
                if outputFile is not None:
                    if isinstance(outputFile, list):
                        self.options['external_output_files'] = outputFile
                    else:
                        self.options['external_output_files'] = [outputFile]
               
                if isinstance(executeCommand, list):
                    self.options['command'] = executeCommand
                else:
                    self.options['command'] = [executeCommand]
                    
                for i in sensitivity.keys():
                    self.deriv_options[i] = sensitivity[i]
                    
                self.changeDir = changeDir
                
                self.saveIteration = saveIteration     
                
                if self.saveIteration:
                    initialInstance(self.analysis.analysisDir)           
                     
            def solve_nonlinear(self, params, unknowns, resids):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysis, self.validInput, self.validGeometry)
               
                # Run parents 
                executeLinked(self.analysis)
                
                if self.saveIteration:
                    storeIteration(self.analysis.analysisDir)
                        
                self.analysis.preAnalysis()
                    
                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise 
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysis.analysisDir)
                    
                #Parent solve_nonlinear function actually runs the external code
                super(openMDAOExternal, self).solve_nonlinear(params, unknowns, resids)
                
                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None: 
                    os.chdir(currentDir)
                    
                self.analysis.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysis.getAnalysisOutVal(varname = i)
            
            # TO add sensitivities       
            #def linearize(self, params, unknowns, resids):
            #    J = {}
            #    return J
        
        class openMDAOComponent(Component):
                
            def __init__(self, 
                         analysis, 
                         inputParam, outputParam, 
                         sensitivity = {}, 
                         saveIteration = False):
                    
                super(openMDAOComponent, self).__init__()
    
                self.analysis = analysis
                    
                # Get check to make sure         
                self.validInput  = self.analysis.getAnalysisVal()
                validOutput      = self.analysis.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysis.capsProblem.geometry.getGeometryVal()
                
                for i in inputParam:
                    
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
                    
                    self.add_param(i, val=value)
                
                for i in outputParam:
                    
                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                   
                for i in sensitivity.keys():
                    self.deriv_options[i] = sensitivity[i]
                    
                self.saveIteration = saveIteration
                
                if self.saveIteration:
                    initialInstance(self.analysis.analysisDir)  
        
            def solve_nonlinear(self, params, unknowns, resids):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysis, self.validInput, self.validGeometry)
               
                if self.saveIteration:
                    storeIteration(self.analysis.analysisDir)
                    
                self.analysis.preAnalysis()
                self.analysis.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysis.output[i].value
            
            # TO add sensitivities       
            #def linearize(self, params, unknowns, resids):
            #    J = {}
            #    return J
    
    #if checkParentExecution(analysis) == False: #Parents can't all self execute
    #    return cCAPS.CAPS_BADVALUE 

    infoDict = analysis.getAnalysisInfo(printInfo = False, infoDict = True)
    
    if executeCommand is not None and infoDict["executionFlag"] == True:

        if versionMajor >=2:
            print("An execution command was provided, but the AIM says it can run itself! Switching ExternalCodeComp to ExplicitComponent")
        else:
            print("An execution command was provided, but the AIM says it can run itself! Switching ExternalCode to Component")
    
    if executeCommand is not None and infoDict["executionFlag"] == False:
        
        if versionMajor >=2:
            return openMDAOExternal_V2(analysis, inputParam, outputParam, 
                                    executeCommand, inputFile, outputFile, sensitivity, changeDir, 
                                    saveIteration)
    
        elif versionMajor == 1 and versionMinor == 7:
            return openMDAOExternal(analysis, inputParam, outputParam, 
                                    executeCommand, inputFile, outputFile, sensitivity, changeDir, 
                                    saveIteration)
    else:
        
        if versionMajor >=2:
           return openMDAOComponent_V2(analysis, inputParam, outputParam, sensitivity, saveIteration)
    
        elif versionMajor == 1 and versionMinor == 7:
            return openMDAOComponent(analysis, inputParam, outputParam, sensitivity, saveIteration)

