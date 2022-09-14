import numpy as np
import math
from corsairlite.core.dataTypes import  TypeCheckedList, Variable
from corsairlite.analysis.models.model import Model

try:
    import torch
except:
    raise ImportError('PyTorch is not provided in the CAPS distribution, and is required for this model.  It is available at https://pytorch.org')


class computeThickness(Model):
    def __init__(self, N=4, inputMode = 0, outputMode = 0, evalLocation = 0.35):
        # Set up all the attributes by calling Model.__init__
        super(computeThickness, self).__init__()
        
        self.description = ('Computes the thickness of an airfoil given a set of Kulfan coefficients.  '+
                            'Various operational modes and transformaions are available by changing self.inputMode.  Order is set by '+
                            'self.N')
        
        self._N = N
        self.evalLocation = evalLocation

        self.inputMode = inputMode
        # Operation Mode 0:  Standard Kulfan
        # Operation Mode 1:  Symmetric Standard Kulfan
        # Operation Mode 2:  Kulfan Transformed by (Au_i + 1), -(Al_i + 1)
        # Operation Mode 3:  Symmetric Kulfan Transformed by (A_i + 1)
        
        self.outputMode = outputMode
        # Mode 0: max thickness
        # Mode 1: evaluated at fixed location
        

        self.outputs.append(Variable("tau" , 0.12,  "", "airfoil thickness to chord ratio"))
        self.availableDerivative = 1
        # self.setupDataContainer()
        
    @property
    def N(self):
        return self._N
    @N.setter
    def N(self,val):
        self._N = val
        self.inputMode = self._inputMode

    @property
    def inputMode(self):
        return self._inputMode
    @inputMode.setter
    def inputMode(self,val):
        self._inputMode = val
        self.inputs = TypeCheckedList(Variable, [])
        if val == 0:
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_%d'%(i+1),  0.2,  '',  'Upper surface Kulfan coefficient, mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('Alower_%d'%(i+1),  -0.2, '',  'Lower surface Kulfan coefficient, mode %d'%(i+1)))
        elif val == 1:
            for i in range(0,self.N):
                self.inputs.append(Variable('A_%d'%(i+1),  0.2,  '',  'Kulfan coefficient, mode %d'%(i+1)))
        elif val == 2:
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_p1_%d'%(i+1),      1.2, '',  'Upper surface Kulfan coefficient, transformed via (A_i + 1), mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('neg_Alower_m1_%d'%(i+1),  1.2, '',  'Lower surface Kulfan coefficient, transformed via -(A_i + 1), mode %d'%(i+1)))
        elif val == 3:
            for i in range(0,self.N):
                self.inputs.append(Variable('A_p1_%d'%(i+1),  1.2, '',  'Kulfan coefficient, transformed via (A_i + 1), mode %d'%(i+1)))
        else:
            raise ValueError('Invalid inputMode')  

    def BlackBox(*args, **kwargs):
        args = list(args)
        self = args.pop(0)
  
        if self.inputMode in [0,2] and len(args) == 2*self.N:
            ipts = args
        elif self.inputMode in [0,2] and len(kwargs.keys()) == 2*self.N:
            ipts = []
            kys = [vr.name for vr in self.inputs]
            for ky in kys:
                ipts.append(kwargs[ky])
        elif self.inputMode in [1,3] and len(args) == self.N:
            ipts = args
        elif self.inputMode in [1,3] and len(kwargs.keys()) == self.N:
            ipts = []
            kys = [vr.name for vr in self.inputs]
            for ky in kys:
                ipts.append(kwargs[ky])
        else:
            raise ValueError('Invalid inputs to the thickness function')
        
        sips = self.sanitizeInputs(*ipts)
        coeffs = np.array([vl.to('').magnitude for vl in sips])
        # These have to stay exactly as inputs for the auto diff to work
        # Make transformation below

        def subfcn(coeffs):
            # Make transformations here
            if self.inputMode == 0:
                Aupper = coeffs[0:self.N]
                Alower = coeffs[self.N:]   
            elif self.inputMode == 1:
                Aupper = coeffs
                Alower = -1 * coeffs    
            elif self.inputMode == 2:
                Aupper = coeffs[0:self.N] -1
                Alower = -1 * coeffs[self.N:] +1
            elif self.inputMode == 3:
                Aupper = coeffs -1
                Alower = -1 * coeffs +1  
            else:
                raise ValueError('Invalid inputMode')     
            
            if self.outputMode == 0:
                ang = np.linspace(0,np.pi/2,101)
                psi = -np.cos(ang)+1
            if self.outputMode == 1:
                psi = np.array([self.evalLocation])
                
            n = len(Aupper) - 1
            N1 = 0.5
            N2 = 1.0

            S = np.zeros([n+1, len(psi)])
            for i in range(0,n+1):
                S[i,:]=math.factorial(n)/(math.factorial(i)*math.factorial(n-i))  *  psi**i  *  (1-psi)**(n-i)

            S = torch.from_numpy(S.astype(float))
            S.requires_grad_(True)

            S_upper = torch.from_numpy(np.zeros(len(psi)).astype(float))
            S_upper.requires_grad_(True)

            S_lower = torch.from_numpy(np.zeros(len(psi)).astype(float))
            S_lower.requires_grad_(True)
            for i in range(0,n+1):
                S_upper = S_upper + Aupper[i] * S[i,:]
                S_lower = S_lower + Alower[i] * S[i,:]

            psi = torch.from_numpy(psi.astype(float))
            psi.requires_grad_(True)
            zeta_lower = psi**(N1)*(1-psi)**(N2)*S_lower
            zeta_upper = psi**(N1)*(1-psi)**(N2)*S_upper

            dist = zeta_upper - zeta_lower
            return max(dist)

        coeffs = torch.from_numpy(coeffs.astype(float))
        coeffs.requires_grad_(True)

        tau = subfcn(coeffs)
        jac = torch.autograd.grad(tau, coeffs, create_graph=True)[0]
        jac = np.array(jac.tolist())

        return float(tau), jac