import os
import copy
import numpy as np
from corsairlite import units
from corsairlite.optimization import Variable
from corsairlite.analysis.models.model import Model
import warnings
            
try:
    import scipy.optimize as spo
except:
    raise ImportError('scipy is not provided in the CAPS distribution, and is required for this model.  It is available at https://scipy.org')

class dragCoefficient(Model):
    def __init__(self):
        super(dragCoefficient, self).__init__()
        
        self.description = 'This model computes the profile drag coefficient based on Hoburg, 2014'
        self.inputs.append( Variable("C_L" ,   1.0,  "", "Wing Lift Coefficient"))
        self.inputs.append( Variable("Re"  ,   1e6,  "", "Reynolds Number"))
        self.inputs.append( Variable("tau" ,  0.12,  "", "Airfoil Thickness to Chord Ratio"))
        self.outputs.append(Variable("C_Dp" , 0.50,  "", "Wing Profile Drag Coefficient"))
        self.availableDerivative = 1
        
        self.setupDataContainer()
    
    def BlackBox(*args, **kwargs):
        args = list(args)       # convert tuple to list
        self = args.pop(0)      # pop off the self argument
        runCases, returnMode, extra = self.parseInputs(*args, **kwargs)
        
        y = []
        dydx = []
        for rc in runCases:
            C_L = rc['C_L'].to('').magnitude
            Re  = rc['Re'].to('').magnitude
            tau = rc['tau'].to('').magnitude

            def residual(x, *args):
                C_Dp = np.exp(x[0])
                C_L = args[0]
                Re = args[1]
                tau = args[2]

                res = (2.56    * C_L**( 5.88) * tau**(-3.32) * Re**(-1.54) * C_Dp**(-2.26) + 
                       3.80e-9 * C_L**(-0.92) * tau**( 6.23) * Re**(-1.38) * C_Dp**(-9.57) + 
                       2.20e-3 * C_L**(-0.01) * tau**( 0.03) * Re**( 0.14) * C_Dp**(-0.73) + 
                       1.19e4  * C_L**( 9.78) * tau**( 1.76) * Re**(-1.00) * C_Dp**(-0.91) + 
                       6.14e-6 * C_L**( 6.53) * tau**(-0.52) * Re**(-0.99) * C_Dp**(-5.19)  -  1)**2            

                return res
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res_implicit = spo.minimize(residual, [np.log(0.05)], args=(C_L,Re,tau), method='SLSQP' )
            
            C_Dp = np.exp(res_implicit['x'][0])
            
            A = [ [2.56    ,  5.88, -3.32, -1.54, -2.26],
                  [3.80e-9 , -0.92,  6.23, -1.38, -9.57],
                  [2.20e-3 , -0.01,  0.03,  0.14, -0.73],
                  [1.19e4  ,  9.78,  1.76, -1.00, -0.91],
                  [6.14e-6 ,  6.53, -0.52, -0.99, -5.19],  ]
            
            Nterms = len(A)
            num_CL = 0.0
            den_CL = 0.0
            num_Re = 0.0
            den_Re = 0.0
            num_tau = 0.0
            den_tau = 0.0
            for i in range(0,Nterms):
                num_CL  += A[i][1] * A[i][0] * C_L**(A[i][1]-1) * tau**A[i][2] * Re**A[i][3] * C_Dp**A[i][4]
                den_CL  += A[i][4] * A[i][0] * C_L**A[i][1]     * tau**A[i][2] * Re**A[i][3] * C_Dp**(A[i][4]-1.0)
                
                num_Re  += A[i][3] * A[i][0] * C_L**A[i][1] * tau**A[i][2] * Re**(A[i][3]-1) * C_Dp**A[i][4]
                den_Re  += A[i][4] * A[i][0] * C_L**A[i][1] * tau**A[i][2] * Re**A[i][3]     * C_Dp**(A[i][4]-1.0)
                
                num_tau += A[i][2] * A[i][0] * C_L**A[i][1] * tau**(A[i][2]-1) * Re**A[i][3] * C_Dp**A[i][4]
                den_tau += A[i][4] * A[i][0] * C_L**A[i][1] * tau**A[i][2]     * Re**A[i][3] * C_Dp**(A[i][4]-1.0)
                
            dCD_dCL  = -1 * num_CL  / den_CL
            dCD_dRe  = -1 * num_Re  / den_Re
            dCD_dtau = -1 * num_tau / den_tau
            
            self.data.addData(dict(zip(['C_L','Re','tau','C_Dp','__D1'],
                                       [C_L, Re, tau, C_Dp, [[dCD_dCL, dCD_dRe, dCD_dtau]]] )))
                              
                
            y.append([C_Dp*units.dimensionless])
            dydx.append([[dCD_dCL*units.dimensionless, dCD_dRe*units.dimensionless, dCD_dtau*units.dimensionless]])
        
        if returnMode < 0:
            returnMode = -1*(returnMode + 1)
            if returnMode == 0:  return y[0][0]
            if returnMode == 1:  return y[0][0], dydx[0][0]
        else:
            if returnMode == 0:
                opt = []
                for i in range(0,len(y)):
                    opt.append( y[i] )
                return opt
            if returnMode == 1:
                opt = []
                for i in range(0,len(y)):
                    opt.append([ y[i],  [dydx[i]] ])
                return opt