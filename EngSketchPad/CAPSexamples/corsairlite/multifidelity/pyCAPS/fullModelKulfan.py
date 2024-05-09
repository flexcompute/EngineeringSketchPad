from corsairlite.core.dataTypes import  TypeCheckedList, Variable
from corsairlite.analysis.wrappers.CAPS.MSES.Kulfan import msesRunKulfan_Alpha
from corsairlite.analysis.models.model import Model
from corsairlite import Q_
import numpy as np

class fullModel(Model):
    def __init__(self, N = 4, inputMode = 0, outputMode = 0, Alpha = None, Re = None, Mach = None,
                         Acrit = None, Airfoil_Points = None, GridAlpha = None, 
                         Coarse_Iteration = None, Fine_Iteration = None,
                         xTransition_Upper = None, xTransition_Lower = None,
                         xGridRange = None, yGridRange = None, cleanup = True, problemObj = None ):
        
        super(fullModel, self).__init__()
        
        self.description = ('This model computes the drag coefficient, lift coefficient, and moment coefficient for a Kulfan defined airfoil'+
                            ' using the CAPS wrapper for MSES.  Various transformations are available with self.inputMode and self.outputMode.  '+
                            ' Order is set by self.N. ')

        self._N = N
        self.inputMode = inputMode
        # Mode 0 : Clean inputs
        # Mode 1 : Clean inputs, no Mach
        # Mode 2 : Clean inputs, geometry only
        # Mode 3 : Transformed inputs
        # Mode 4 : Transformed inputs, no Mach
        # Mode 5 : Transformed inputs, geometry only
        
        self.outputMode = outputMode
        # Mode 0 : Clean outputs
        # Mode 1 : Clean CL, CD
        # Mode 2 : Drag Coefficient only
        # Mode 3 : Clean CD, CL, transformed CM
        
        self.Alpha = Alpha
        self.Re = Re
        self.Mach = Mach

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
        self.problemObj = problemObj


    @property
    def N(self):
        return self._N
    @N.setter
    def N(self,val):
        self._N = val
        self.inputMode = self._inputMode
        self.outputMode = self._outputMode

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
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_%d'%(i+1),  0.2, '',  'Upper surface Kulfan coefficient, mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('Alower_%d'%(i+1), -0.2, '',  'Lower surface Kulfan coefficient, mode %d'%(i+1)))
        elif self.inputMode == 1:
            self.inputs.append( Variable('Alpha',          2.0,  'deg', 'Angle of Attack'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_%d'%(i+1),  0.2, '',  'Upper surface Kulfan coefficient, mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('Alower_%d'%(i+1), -0.2, '',  'Lower surface Kulfan coefficient, mode %d'%(i+1)))
        elif self.inputMode == 2:
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_%d'%(i+1),  0.2, '',  'Upper surface Kulfan coefficient, mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('Alower_%d'%(i+1), -0.2, '',  'Lower surface Kulfan coefficient, mode %d'%(i+1)))
        elif self.inputMode == 3:
            self.inputs.append( Variable('Alpha',          22.0,  'deg', 'Angle of Attack + 20 deg'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            self.inputs.append( Variable('Mach',           0.2,      '', 'Mach Number'))
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_p1_%d'%(i+1),      1.2, '',  'Upper surface Kulfan coefficient, transformed via (A_i + 1), mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('neg_Alower_m1_%d'%(i+1),  1.2, '',  'Lower surface Kulfan coefficient, transformed via -(A_i + 1), mode %d'%(i+1)))
        elif self.inputMode == 4:
            self.inputs.append( Variable('Alpha',          22.0,  'deg', 'Angle of Attack + 20 deg'))
            self.inputs.append( Variable('Re',             1e6,      '', 'Reynolds Number'))
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_p1_%d'%(i+1),      1.2, '',  'Upper surface Kulfan coefficient, transformed via (A_i + 1), mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('neg_Alower_m1_%d'%(i+1),  1.2, '',  'Lower surface Kulfan coefficient, transformed via -(A_i + 1), mode %d'%(i+1)))
        elif self.inputMode == 5:
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_p1_%d'%(i+1),      1.2, '',  'Upper surface Kulfan coefficient, transformed via (A_i + 1), mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('neg_Alower_m1_%d'%(i+1),  1.2, '',  'Lower surface Kulfan coefficient, transformed via -(A_i + 1), mode %d'%(i+1)))
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
    
        if self.inputMode in [0,1,3,4]:
            Alpha = sips.pop(0).to('deg').magnitude
            if self.inputMode in [3,4]:
                Alpha -= 20.0
        else:
            if self.Alpha is not None:
                Alpha = self.Alpha
            else:
                raise ValueError('Must set self.Alpha in this inputMode')     
            
        if self.inputMode in [0,1,3,4]:
            Re = sips.pop(0).to('').magnitude
        else:
            if self.Re is not None:
                Re = self.Re
            else:
                raise ValueError('Must set self.Re in this inputMode')       
            
        if self.inputMode in [0,3]:
            Mach = sips.pop(0).to('').magnitude
        else:
            if self.Mach is not None:
                Mach = self.Mach
            else:
                raise ValueError('Must set self.Mach in this inputMode')
                
        Avals = [vl.to('').magnitude for vl in sips]
        Aupper = np.array(Avals[0:self.N])
        Alower = np.array(Avals[self.N:])
        
        if self.inputMode in [3,4,5]:
            Aupper =     Aupper - 1
            Alower = -1*(Alower - 1)
            
        opt = msesRunKulfan_Alpha(Aupper = Aupper, Alower = Alower,
                                 Mach = Mach, Re = Re, Alpha = Alpha,
                                 Acrit = self.Acrit, Airfoil_Points = self.Airfoil_Points, GridAlpha = self.GridAlpha, 
                                 Coarse_Iteration = self.Coarse_Iteration, Fine_Iteration = self.Fine_Iteration,
                                 xTransition_Upper = self.xTransition_Upper, xTransition_Lower = self.xTransition_Lower,
                                 xGridRange = self.xGridRange, yGridRange = self.yGridRange, cleanup = self.cleanup,
                                 problemObj = self.problemObj )

        CD      =    opt[0]['CD']
        CL      =    opt[0]['CL']
        CM      =    opt[0]['CM']

        CD_Alpha     = opt[1]['CD']['Alpha']
        CD_Re        = opt[1]['CD']['Re']
        CD_Mach      = opt[1]['CD']['Mach']
        CD_aupper    = opt[1]['CD']['aupper']
        CD_alower    = opt[1]['CD']['alower']
        
        CL_Alpha     = opt[1]['CL']['Alpha']
        CL_Re        = opt[1]['CL']['Re']
        CL_Mach      = opt[1]['CL']['Mach']
        CL_aupper    = opt[1]['CL']['aupper']
        CL_alower    = opt[1]['CL']['alower']
        
        CM_Alpha     = opt[1]['CM']['Alpha']
        CM_Re        = opt[1]['CM']['Re']
        CM_Mach      = opt[1]['CM']['Mach']
        CM_aupper    = opt[1]['CM']['aupper']
        CM_alower    = opt[1]['CM']['alower']
        
        if self.inputMode in [3,4,5]:
            CD_alower = [vl*-1 for vl in CD_alower]
            CL_alower = [vl*-1 for vl in CL_alower]
            CM_alower = [vl*-1 for vl in CM_alower]   
            
        if self.outputMode == 3:
            CM_Alpha     *= -1
            CM_Re        *= -1
            CM_Mach      *= -1
            CM_aupper    = [vl*-1 for vl in CM_aupper]
            CM_alower    = [vl*-1 for vl in CM_alower]
        
        if self.inputMode in [0,3]:
            CD_ders = [CD_Alpha , CD_Re, CD_Mach] + CD_aupper + CD_alower
            CL_ders = [CL_Alpha , CL_Re, CL_Mach] + CL_aupper + CL_alower
            CM_ders = [CM_Alpha , CM_Re, CM_Mach] + CM_aupper + CM_alower
        elif self.inputMode in [1,4]:
            CD_ders = [CD_Alpha , CD_Re] + CD_aupper + CD_alower
            CL_ders = [CL_Alpha , CL_Re] + CL_aupper + CL_alower
            CM_ders = [CM_Alpha , CM_Re] + CM_aupper + CM_alower
        elif self.inputMode in [2,5]:
            CD_ders = CD_aupper + CD_alower
            CL_ders = CL_aupper + CL_alower
            CM_ders = CM_aupper + CM_alower
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
