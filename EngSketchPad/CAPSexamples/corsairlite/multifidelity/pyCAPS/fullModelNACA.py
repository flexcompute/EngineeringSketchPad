from corsairlite.core.dataTypes import  TypeCheckedList, Variable
from corsairlite.analysis.wrappers.CAPS.MSES.NACA4 import msesRunNACA4_Alpha
from corsairlite.analysis.models.model import Model
from corsairlite import Q_

class fullModel(Model):
    def __init__(self, inputMode = 0, outputMode = 0, Alpha = None, Re = None, Mach = None,
                         maxCamber = None, camberDistance = None,
                         Acrit = None, Airfoil_Points = None, GridAlpha = None, 
                         Coarse_Iteration = None, Fine_Iteration = None,
                         xTransition_Upper = None, xTransition_Lower = None,
                         xGridRange = None, yGridRange = None, cleanup = True, problemObj = None ):
        super(fullModel, self).__init__()
        
        self.description = ('This model computes the drag coefficient, lift coefficient, and moment coefficient for a NACA 4 series airfoil'+
                            ' using the CAPS wrapper for MSES.  Various transformations are available with self.inputMode and self.outputMode.')
        

        self.inputMode = inputMode
        # Mode 0 : Clean inputs
        # Mode 1 : Clean inputs, no Mach
        # Mode 2 : Clean inputs alpha, Re, thickness only
        # Mode 3 : Clean inputs, geometry only
        # Mode 4 : Transformed alpha
        # Mode 5 : Transformed alpha, no Mach
        # Mode 6 : Transformed input alpha, clean Re, clean thickness only
        
        self.outputMode = outputMode
        # Mode 0 : Clean outputs
        # Mode 1 : Clean CL, CD
        # Mode 2 : Drag Coefficient only
        # Mode 3 : Clean CD, CM, transformed CM

        
        self.Alpha = Alpha
        self.Re = Re
        self.Mach = Mach
        self.maxCamber = maxCamber
        self.camberDistance = camberDistance

        self.Acrit             =  Acrit            
        self.Airfoil_Points    =  Airfoil_Points   
        self.GridAlpha         =  GridAlpha        
        self.Coarse_Iteration  =  Coarse_Iteration 
        self.Fine_Iteration    =  Fine_Iteration   
        self.xTransition_Upper =  xTransition_Upper
        self.xTransition_Lower =  xTransition_Lower
        self.xGridRange        =  xGridRange       
        self.yGridRange        =  yGridRange       
        self.cleanup           =  cleanup          
        
        self.availableDerivative = 1
        # self.setupDataContainer()

        self.problemObj = problemObj

    @property
    def inputMode(self):
        return self._inputMode
    @inputMode.setter
    def inputMode(self,val):
        self._inputMode = val
        self.inputs = TypeCheckedList(Variable, [])
        if self.inputMode == 0:
            self.inputs.append( Variable('Alpha',          2.0,   'deg', 'Angle of Attack'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            self.inputs.append( Variable('Mach',           0.2,      '', 'Mach Number'))
            self.inputs.append( Variable('maxCamber',      0.02,     '', 'Maximum camber scaled by chord, can be any float number'))
            self.inputs.append( Variable('camberDistance', 0.4,      '', 'Distance along chord where max camber occours scaled by chord, can be any float number'))
            self.inputs.append( Variable('thickness',      0.12,     '', 'Maximum thickness of the airfoil scaled by chord, can be any float number'))
        elif self.inputMode == 1:
            self.inputs.append( Variable('Alpha',          2.0,  'deg', 'Angle of Attack'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            self.inputs.append( Variable('maxCamber',      0.02,     '', 'Maximum camber scaled by chord, can be any float number'))
            self.inputs.append( Variable('camberDistance', 0.4,      '', 'Distance along chord where max camber occours scaled by chord, can be any float number'))
            self.inputs.append( Variable('thickness',      0.12,     '', 'Maximum thickness of the airfoil scaled by chord, can be any float number'))
        elif self.inputMode == 2:
            self.inputs.append( Variable('Alpha',          2.0,  'deg', 'Angle of Attack'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            self.inputs.append( Variable('thickness',      0.12,     '', 'Maximum thickness of the airfoil scaled by chord, can be any float number'))
        elif self.inputMode == 3:
            self.inputs.append( Variable('maxCamber',      0.02,     '', 'Maximum camber scaled by chord, can be any float number'))
            self.inputs.append( Variable('camberDistance', 0.4,      '', 'Distance along chord where max camber occours scaled by chord, can be any float number'))
            self.inputs.append( Variable('thickness',      0.12,     '', 'Maximum thickness of the airfoil scaled by chord, can be any float number'))
        elif self.inputMode == 4:
            self.inputs.append( Variable('Alpha',          22.0,  'deg', 'Angle of Attack + 20 deg'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            self.inputs.append( Variable('Mach',           0.2,      '', 'Mach Number'))
            self.inputs.append( Variable('maxCamber',      0.02,     '', 'Maximum camber scaled by chord, can be any float number'))
            self.inputs.append( Variable('camberDistance', 0.4,      '', 'Distance along chord where max camber occours scaled by chord, can be any float number'))
            self.inputs.append( Variable('thickness',      0.12,     '', 'Maximum thickness of the airfoil scaled by chord, can be any float number'))
        elif self.inputMode == 5:
            self.inputs.append( Variable('Alpha',          22.0,  'deg', 'Angle of Attack + 20 deg'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            self.inputs.append( Variable('maxCamber',      0.02,     '', 'Maximum camber scaled by chord, can be any float number'))
            self.inputs.append( Variable('camberDistance', 0.4,      '', 'Distance along chord where max camber occours scaled by chord, can be any float number'))
            self.inputs.append( Variable('thickness',      0.12,     '', 'Maximum thickness of the airfoil scaled by chord, can be any float number'))
        elif self.inputMode == 6:
            self.inputs.append( Variable('Alpha',          22.0,  'deg', 'Angle of Attack + 20 deg'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            self.inputs.append( Variable('thickness',      0.12,     '', 'Maximum thickness of the airfoil scaled by chord, can be any float number'))    
        else:
            raise ValueError('Invalid inputMode')

    @property
    def outputMode(self):
        return self._outputMode
    @outputMode.setter
    def outputMode(self,val):
        self._outputMode = val
        self.outputs = TypeCheckedList(Variable, [])
        if self.outputMode == 0:
            self.outputs.append(Variable('C_D',            0.05,     '', 'Profile Drag Coefficient'))
            self.outputs.append(Variable('C_L',             0.5,     '', 'Lift Coefficient'))
            self.outputs.append(Variable('C_M',           -0.05,     '', 'Moment Coefficient'))
        elif self.outputMode == 1:
            self.outputs.append(Variable('C_D',            0.05,     '', 'Profile Drag Coefficient'))
            self.outputs.append(Variable('C_L',            0.5,      '', 'Lift Coefficient'))
        elif self.outputMode == 2:
            self.outputs.append(Variable('C_D',            0.05,     '', 'Profile Drag Coefficient'))
        elif self.outputMode == 3:
            self.outputs.append(Variable('C_D',            0.05,     '', 'Profile Drag Coefficient'))
            self.outputs.append(Variable('C_L',            0.5,      '', 'Lift Coefficient'))
            self.outputs.append(Variable('negCMp1',        1.05,     '', 'Negative moment coefficient Plus 1'))
        else:
            raise ValueError('Invalid outputMode')
            
    @property
    def Alpha(self):
        return self._Alpha
    @Alpha.setter
    def Alpha(self,val):
        if isinstance(val, Q_):
            self._Alpha = val.to('deg').magnitude
        else:
            self._Alpha = val

    @property
    def Re(self):
        return self._Re
    @Re.setter
    def Re(self,val):
        if isinstance(val, Q_):
            self._Re = val.to('').magnitude
        else:
            self._Re = val
            
    @property
    def Mach(self):
        return self._Mach
    @Mach.setter
    def Mach(self,val):
        if isinstance(val, Q_):
            self._Mach = val.to('').magnitude
        else:
            self._Mach = val
            
    @property
    def maxCamber(self):
        return self._maxCamber
    @maxCamber.setter
    def maxCamber(self,val):
        if isinstance(val, Q_):
            self._maxCamber = val.to('').magnitude
        else:
            self._maxCamber = val

    @property
    def camberDistance(self):
        return self._camberDistance
    @camberDistance.setter
    def camberDistance(self,val):
        if isinstance(val, Q_):
            self._camberDistance = val.to('').magnitude
        else:
            self._camberDistance = val
            

    def BlackBox(*args, **kwargs):
        args = list(args)
        self = args.pop(0)
        
        if len(args) + len(kwargs.keys()) != len(self.inputs):
            raise ValueError('Incorrect number of inputs')
        
        kys = [vr.name for vr in self.inputs]
        trimKys = kys[0:len(args)]
        argDict = dict(zip(trimKys, args))
        iptsDict = argDict | kwargs
        ipts = [iptsDict[nm] for nm in kys]
    
        sips = self.sanitizeInputs(*ipts)
            
        if self.inputMode in [0,1,2,4,5,6]:
            Alpha = sips.pop(0).to('deg').magnitude
            if self.inputMode in [4,5,6]:
                Alpha -= 20.0
        else:
            if self.Alpha is not None:
                Alpha = self.Alpha
            else:
                raise ValueError('Must set self.Alpha in this inputMode')     
            
        if self.inputMode in [0,1,2,4,5,6]:
            Re = sips.pop(0).to('').magnitude
        else:
            if self.Re is not None:
                Re = self.Re
            else:
                raise ValueError('Must set self.Re in this inputMode')       
            
        if self.inputMode in [0,4]:
            Mach = sips.pop(0).to('').magnitude
        else:
            if self.Mach is not None:
                Mach = self.Mach
            else:
                raise ValueError('Must set self.Mach in this inputMode')
                
        if self.inputMode in [0,1,3,4,5]:
            maxCamber = sips.pop(0).to('').magnitude
        else:
            if self.maxCamber is not None:
                maxCamber = self.maxCamber
            else:
                raise ValueError('Must set self.maxCamber in this inputMode')
                
        if self.inputMode in [0,1,3,4,5]:
            camberDistance = sips.pop(0).to('').magnitude
        else:
            if self.camberDistance is not None:
                camberDistance = self.camberDistance
            else:
                raise ValueError('Must set self.camberDistance in this inputMode')
                
        thickness = sips[0].to('').magnitude
                              
        opt = msesRunNACA4_Alpha(camber=maxCamber,  maxloc=camberDistance, thickness=thickness,
                                 Mach = Mach, Re = Re, Alpha = Alpha,
                                 Acrit = self.Acrit, Airfoil_Points = self.Airfoil_Points, GridAlpha = self.GridAlpha, 
                                 Coarse_Iteration = self.Coarse_Iteration, Fine_Iteration = self.Fine_Iteration,
                                 xTransition_Upper = self.xTransition_Upper, xTransition_Lower = self.xTransition_Lower,
                                 xGridRange = self.xGridRange, yGridRange = self.yGridRange, cleanup = self.cleanup,
                                 problemObj = self.problemObj )
                            
            
        CD      =    opt[0]['CD']
        CL      =    opt[0]['CL']
        CM      =    opt[0]['CM']

        CD_Alpha          = opt[1]['CD']['Alpha']
        CD_Re             = opt[1]['CD']['Re']
        CD_Mach           = opt[1]['CD']['Mach']
        CD_maxCamber      = opt[1]['CD']['camber']
        CD_camberDistance = opt[1]['CD']['maxloc']
        CD_thickness      = opt[1]['CD']['thickness']
        
        CL_Alpha          = opt[1]['CL']['Alpha']
        CL_Re             = opt[1]['CL']['Re']
        CL_Mach           = opt[1]['CL']['Mach']
        CL_maxCamber      = opt[1]['CL']['camber']
        CL_camberDistance = opt[1]['CL']['maxloc']
        CL_thickness      = opt[1]['CL']['thickness']
        
        CM_Alpha          = opt[1]['CM']['Alpha']
        CM_Re             = opt[1]['CM']['Re']
        CM_Mach           = opt[1]['CM']['Mach']
        CM_maxCamber      = opt[1]['CM']['camber']
        CM_camberDistance = opt[1]['CM']['maxloc']
        CM_thickness      = opt[1]['CM']['thickness']
        
        if self.outputMode == 3:
            CM_Alpha          *= -1
            CM_Re             *= -1
            CM_Mach           *= -1
            CM_maxCamber      *= -1
            CM_camberDistance *= -1
            CM_thickness      *= -1
        
        if self.inputMode in [0,4]:
            CD_ders = [CD_Alpha , CD_Re, CD_Mach, CD_maxCamber, CD_camberDistance, CD_thickness]
            CL_ders = [CL_Alpha , CL_Re, CL_Mach, CL_maxCamber, CL_camberDistance, CL_thickness]
            CM_ders = [CM_Alpha , CM_Re, CM_Mach, CM_maxCamber, CM_camberDistance, CM_thickness]
        elif self.inputMode in [1,5]:
            CD_ders = [CD_Alpha , CD_Re, CD_maxCamber, CD_camberDistance, CD_thickness]
            CL_ders = [CL_Alpha , CL_Re, CL_maxCamber, CL_camberDistance, CL_thickness]
            CM_ders = [CM_Alpha , CM_Re, CM_maxCamber, CM_camberDistance, CM_thickness]
        elif self.inputMode in [2,6]:
            CD_ders = [CD_Alpha , CD_Re, CD_thickness]
            CL_ders = [CL_Alpha , CL_Re, CL_thickness]
            CM_ders = [CM_Alpha , CM_Re, CM_thickness]
        elif self.inputMode in [3]:
            CD_ders = [CD_maxCamber, CD_camberDistance, CD_thickness]
            CL_ders = [CL_maxCamber, CL_camberDistance, CL_thickness]
            CM_ders = [CM_maxCamber, CM_camberDistance, CM_thickness]
        else:
            raise ValueError('Invalid inputMode')
        
        if self.outputMode == 0:
            outputVals = [CD, CL, CM]
            outputDers = [CD_ders, CL_ders, CM_ders]
        elif self.outputMode == 1:
            outputVals = [CD, CL]
            outputDers = [CD_ders, CL_ders]
        elif self.outputMode == 2:
            outputVals = CD
            outputDers = CD_ders
        elif self.outputMode == 3:
            negCMp1 = -1*CM + 1
            outputVals = [CD, CL, negCMp1]
            outputDers = [CD_ders, CL_ders, CM_ders]
        else:
            raise ValueError('Invalid outputMode')
        return outputVals, outputDers
