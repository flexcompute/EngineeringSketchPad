import numpy as np
import math
from corsairlite.core.dataTypes.variable import Variable
from corsairlite.core.dataTypes.kulfan import Kulfan
from corsairlite.analysis.models.model import Model

try:
    import torch
except:
    raise ImportError('PyTorch is not provided in the CAPS distribution, and is required for this model.  It is available at https://pytorch.org')


def cby(IFUN_in,S,Nmodes):
    if S < 0.0:
        sgnFlp = -1.0
        S = abs(S)
    else:
        sgnFlp = 1.0
        
    if IFUN_in > Nmodes/2:
        IFUN = IFUN_in-Nmodes/2
    else:
        IFUN = IFUN_in
        
    SWT = 2.0
    SW = SWT*S / (1.0 + (SWT-1)*S)
    X = 1.0 - 2.0*SW
    X = min([1.0,X])
    X = max([-1.0,X])
    THETA = np.arccos(X)
    RF = float(IFUN + 1)
    if IFUN % 2 == 0:
        GFUNTP = (  X - np.cos(RF*THETA))/RF
    else:
        GFUNTP = (1.0 - np.cos(RF*THETA))/RF
        
    fullMode = sgnFlp*GFUNTP
    # fullMode = GFUNTP
    if IFUN_in > Nmodes/2:
        if sgnFlp>0:
            fullMode = 0.0   
    else:
        if sgnFlp < 0:
            fullMode = 0.0
    return fullMode

def computeCoordinates(wts, baselineAirfoil=None):
    if baselineAirfoil is None:
        baselineAirfoil = Kulfan()
        baselineAirfoil.naca4_like(0,0,12)

    xcoords = baselineAirfoil.xcoordinates
    ycoords = baselineAirfoil.ycoordinates

    sArcs = []
    sArc = 0.0
    LEix = None
    for i in range(0,len(xcoords)):
        if i==0:
            sArcs.append(0.0)
        else:
            dxl = xcoords[i] - xcoords[i-1]
            dyl = ycoords[i] - ycoords[i-1]
            if xcoords[i] == 0.0 and ycoords[i] == 0.0:
                LEix = i
            sArc += ((dxl)**2 + (dyl)**2)**0.5
            sArcs.append(sArc)

    upperArcs = np.array(sArcs[0:LEix+1])
    lowerArcs = np.array(sArcs[LEix:])

    upperArcs = 1-upperArcs/upperArcs[-1]
    lowerArcs = -1 * (lowerArcs-lowerArcs[0])/(lowerArcs[-1]-lowerArcs[0])
    sArcs = np.append(upperArcs[0:-1],lowerArcs)

    import torch
    N = len(wts)
    Nmodes = N
    modes = np.linspace(1,N,N)
    polys = []
    for md in modes:
        vls = []
        svls = sArcs
        for s in svls:
            vls.append(cby(md,s,Nmodes))
        polys.append(np.array(vls))
    scaledPolys = torch.zeros([N, len(svls)])
    for i in range(0,N):
        scaledPolys[i,:] = wts[i] * torch.from_numpy(polys[i])
    cheby_vals = torch.sum(scaledPolys, 0)

    movedX  = []
    movedY  = []

    for i in range(0,len(xcoords)):
        if i==0:
            movedX.append(xcoords[i])
            movedY.append(ycoords[i])
        elif i != len(xcoords)-1:
            dxl = xcoords[i] - xcoords[i-1]
            dxr = xcoords[i+1] - xcoords[i]
            dyl = ycoords[i] - ycoords[i-1]
            dyr = ycoords[i+1] - ycoords[i]

            norm_app = np.array([dyr + dyl , -1 * (dxr + dxl) ])
            norm_app = norm_app/(norm_app[0]**2 + norm_app[1]**2)**0.5

            if sArcs[i] >= 0 :
                movedX.append( xcoords[i]+cheby_vals[i]*norm_app[0])
                movedY.append( ycoords[i]+cheby_vals[i]*norm_app[1])
            else:
                movedX.append( xcoords[i]-cheby_vals[i]*norm_app[0])
                movedY.append( ycoords[i]-cheby_vals[i]*norm_app[1])
        else:
            movedX.append(xcoords[i])
            movedY.append(ycoords[i])

    return [movedX, movedY]

class computeThickness(Model):
    def __init__(self, N = 20, inputMode = 1, baselineAirfoil = None):
        # Set up all the attributes by calling Model.__init__
        super(computeThickness, self).__init__()
        
        self.description = ('Computes the thickness of an airfoil from the MSES Chebyshev modes.  '+ 
                            'Various transformations are available by setting self.inputMode.  ' + 
                            'Number of modes is set via input N=[Nmodes] or via self.N=Nmodes.  ' + 
                            'Baseline airfoil defaults to 0012, but can be set via self.baselineAirfoil ' +
                            'or passed in as an input to baselineAirfoil')

        self._N = N
        self.inputMode = inputMode
        # Mode 0 : Standard MSES definition 
        # Mode 1 : plus one shift
        self.baselineAirfoil = baselineAirfoil
        
        self.outputs.append(Variable("tau" , 0.12,  "", "airfoil thickness to chord ratio"))
        self.availableDerivative = 1
        self.setupDataContainer()

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
                self.inputs.append(Variable('Aupper_%d'%(i+1),   0.02, '',  'Upper surface Chebyshev weight, mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('Alower_%d'%(i+1),   0.02, '',  'Lower surface Chebyshev weight, mode %d'%(i+1)))
        elif val == 1:
            for i in range(0,self.N):
                self.inputs.append(Variable('Aupper_p1_%d'%(i+1),   1.02, '',  'Upper surface Chebyshev weight, transformed via (A_i + 1), mode %d'%(i+1)))
            for i in range(0,self.N):
                self.inputs.append(Variable('Alower_p1_%d'%(i+1),   1.02, '',  'Lower surface Chebyshev weight, transformed via (A_i + 1), mode %d'%(i+1)))
        else:
            raise ValueError('Invalid inputMode')  

    def BlackBox(*args, **kwargs):
        args = list(args)
        self = args.pop(0)
        
        if len(args) == 2*self.N:
            ipts = args
        elif len(kwargs.keys()) == 2*self.N:
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
                wts = coeffs
            elif self.inputMode == 1:
                wts = coeffs-1
            else:
                raise ValueError('Invalid inputMode')

            if self.baselineAirfoil is None:
                baselineAirfoil = Kulfan()
                baselineAirfoil.naca4_like(0,0,12)
            else:
                baselineAirfoil = self.baselineAirfoil

            xcoords = baselineAirfoil.xcoordinates
            ycoords = baselineAirfoil.ycoordinates

            sArcs = []
            sArc = 0.0
            for i in range(0,len(xcoords)):
                if i==0:
                    sArcs.append(0.0)
                else:
                    dxl = xcoords[i] - xcoords[i-1]
                    dyl = ycoords[i] - ycoords[i-1]
                    sArc += ((dxl)**2 + (dyl)**2)**0.5
                    sArcs.append(sArc)
            sArcs = np.array(sArcs)/sArcs[-1] * 2
            sArcs -= 1.0
            sArcs *= -1
            splitIndex = np.where(sArcs<=0)[0][0]
            
            N = len(wts)
            Nmodes = N
            modes = np.linspace(1,N,N)
            polys = []
            for md in modes:
                vls = []
                svls = sArcs
                for s in svls:
                    vls.append(cby(md,s,Nmodes))
                polys.append(np.array(vls))
            scaledPolys = torch.zeros([N, len(svls)])
            
            #start tracking now

            for i in range(0,N):
                scaledPolys[i,:] = wts[i] * torch.from_numpy(polys[i])
            cheby_vals = torch.sum(scaledPolys, 0)

            movedX  = []
            movedY  = []

            for i in range(0,len(xcoords)):
                if i==0:
                    movedX.append(xcoords[i])
                    movedY.append(ycoords[i])
                elif i != len(xcoords)-1:
                    dxl = xcoords[i] - xcoords[i-1]
                    dxr = xcoords[i+1] - xcoords[i]
                    dyl = ycoords[i] - ycoords[i-1]
                    dyr = ycoords[i+1] - ycoords[i]

                    norm_app = np.array([dyr + dyl , -1 * (dxr + dxl) ])
                    norm_app = norm_app/(norm_app[0]**2 + norm_app[1]**2)**0.5

                    if sArcs[i] >= 0 :
                        movedX.append( xcoords[i]+cheby_vals[i]*norm_app[0])
                        movedY.append( ycoords[i]+cheby_vals[i]*norm_app[1])
                    else:
                        movedX.append( xcoords[i]-cheby_vals[i]*norm_app[0])
                        movedY.append( ycoords[i]-cheby_vals[i]*norm_app[1])
                else:
                    movedX.append(xcoords[i])
                    movedY.append(ycoords[i])

            ang = torch.linspace(0,np.pi/2,301)
            psi = -torch.cos(ang)+1
            
            xUpper = list(reversed(movedX[0:splitIndex]))
            xLower = movedX[splitIndex:]

            yUpper = list(reversed(movedY[0:splitIndex]))
            yLower = movedY[splitIndex:]
                
            tauStack = [None] * len(psi)
            for i in range(0,len(psi)):
                psi_coord = psi[i]
                if psi_coord == 0.0 or psi_coord == 1.0:
                    tauStack[i] = 0.0
                else:
                    try:
                        for ii in range(0,len(xUpper)):
                            if xUpper[ii] >= psi_coord:
                                upper_ux = ii
                                break
                        if xUpper[upper_ux] == psi_coord:
                            yUpper_interp = yUpper[upper_ux]
                        else:
                            upper_lx = upper_ux - 1
                            if upper_lx < 0:
                                raise ValueError('Lower Index less than zero')
                            deltaY = yUpper[upper_ux] - yUpper[upper_lx]
                            deltaX = xUpper[upper_ux] - xUpper[upper_lx]
                            slp = deltaY/deltaX
                            subDeltaX = psi_coord - xUpper[upper_lx]
                            yUpper_interp = yUpper[upper_lx] + slp*subDeltaX

                        for ii in range(0,len(xLower)):
                            if xLower[ii] >= psi_coord:
                                lower_ux = ii
                                break
                        if xLower[lower_ux] == psi_coord:
                            yLower_interp = yLower[lower_ux]
                        else:
                            lower_lx = lower_ux - 1
                            if lower_lx < 0:
                                raise ValueError('Lower Index less than zero')
                            deltaY = yLower[lower_ux] - yLower[lower_lx]
                            deltaX = xLower[lower_ux] - xLower[lower_lx]
                            slp = deltaY/deltaX
                            subDeltaX = psi_coord - xLower[lower_lx]
                            yLower_interp = yLower[lower_lx] + slp*subDeltaX
                        
                        tauStack[i] = yUpper_interp - yLower_interp
                    except:
                        tauStack[i] = 0

            return max(tauStack)

        coeffs = torch.from_numpy(coeffs.astype(float))
        coeffs.requires_grad_(True)

        tau = subfcn(coeffs)
        jac = torch.autograd.grad(tau, coeffs, create_graph=True)[0]
        jac = np.array(jac.tolist())
        
        return float(tau), jac