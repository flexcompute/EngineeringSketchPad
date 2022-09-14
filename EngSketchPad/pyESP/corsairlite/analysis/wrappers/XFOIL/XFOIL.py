import subprocess
import numpy as np
from corsairlite.core.dataTypes import Kulfan
from corsairlite import units
import os
import math
import shutil
path_to_XFOIL = shutil.which('xfoil')

def XFOIL(mode, 
          airfoil,
          val = 0.0, 
          minval = 0.0,
          maxval = 1.0,
          steps = 10,
          Re = 1e7,
          # M = 0.0,
          flapLocation = None,
          flapDeflection = 0.0,
          polarfile='polars.txt',
          defaultDatfile = 'airfoil.dat',
          saveDatfile = False,
          removeData = True,
          max_iter=100):

    M = 0.0
    
    mode = mode.lower()
    if mode == 'alpha':
        mode = 'alfa'

    if mode not in ['alfa','cl','aseq','cseq']:
        raise ValueError('Invalid input mode.  Must be one of: alfa, cl, aseq, cseq.')
    
    if os.path.isfile(polarfile):
        os.remove(polarfile)

    if isinstance(airfoil, str):        
        if ('.dat' in airfoil) or ('.txt' in airfoil):
            if os.path.isfile(airfoil):
                topline = 'load ' + airfoil + ' \n afl \n'
            else:
                raise ValueError('Could not find airfoil to be read')
        ck1 = 'naca' == airfoil.lower()[0:4]
        ck2 = airfoil[-4:].isdigit()
        if ck1:
            if ck2 and (len(airfoil)!=8):
                afl = airfoil.split()
                airfoil = afl[0]+afl[1]
            if ck2 and (len(airfoil)==8):
                topline = airfoil + ' \n'
            else:
                raise ValueError('Could not parse the NACA 4 digit airfoil')
        
    elif isinstance(airfoil, Kulfan):
        if os.path.isfile(defaultDatfile):
            os.remove(defaultDatfile)
        airfoil.write2file(defaultDatfile)
        topline = 'load ' + defaultDatfile + ' \n' + 'airfoil \n'
        
    else:
        raise ValueError('Could not parse airfoil')
    
    proc = subprocess.Popen([path_to_XFOIL], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    estr = ''
    estr += topline
    estr += 'plop\n'
    estr += 'g\n'
    estr += '\n'
    if flapLocation is not None:
        ck1 = flapLocation >= 0.0
        ck2 = flapLocation <= 1.0
        if ck1 and ck2:
            estr += 'gdes \n'
            estr += 'flap \n'
            estr += '%f \n'%(flapLocation)
            estr += '999 \n'
            estr += '0.5 \n'
            estr += '%f \n'%(flapDeflection)
            estr += 'x \n'
            estr += '\n'
        else:
            raise ValueError('Invalid flapLocation.  Must be between 0.0 and 1.0')
    estr += 'oper \n'
    estr += "iter %d\n" %(max_iter)
    estr += 'visc \n'
    estr += "%.0f \n" %(float(Re))
    estr += "M \n"
    estr += "%.2f \n" %(M)
    if mode == 'aseq':
        amin = math.floor(minval)
        estr += "aseq %.2f %.2f %.2f\n" %(0.0, amin, 1.0)
    estr += 'pacc \n'
    estr += '\n'
    estr += '\n'
    if mode == 'alfa':
        estr += "alfa %.2f \n" %(val)
    if mode == 'cl':
        estr += "cl %.3f \n" %(val)
    if mode == 'aseq':
        amin = math.floor(minval)
        amax = math.ceil(maxval)
        estr += "aseq %.2f %.2f %.2f\n" %(amin, amax, 1.0)
    if mode == 'cseq':
        estr += "cseq %.3f %.3f, %.3f\n" %(minval, maxval, (maxval-minval)/float(steps))
    estr += 'pwrt \n'
    estr += polarfile + ' \n'
    estr += '\n'
    estr += 'q \n'
    f = open('wrtfl.txt','w')
    f.write(estr)
    f.close()
    proc.stdin.write(estr.encode())
    stdout_val = proc.communicate()[0]
    proc.stdin.close()

    data = np.genfromtxt(polarfile, skip_header=12)

    if mode == 'alfa' or mode == 'cl':
        alpha   = data[0]
        cl      = data[1]
        cd      = data[2]
        cdp     = data[3]
        cm      = data[4]
        xtr_top = data[5]
        xtr_bot = data[6]
        Reval   = Re
        Mval    = M
    if mode == 'aseq' or mode == 'cseq':
        alpha   = data[:,0]
        cl      = data[:,1]
        cd      = data[:,2]
        cdp     = data[:,3]
        cm      = data[:,4]
        xtr_top = data[:,5]
        xtr_bot = data[:,6]
        Reval   = Re * np.ones(len(alpha))
        Mval    = M * np.ones(len(alpha))

    if os.path.isfile(':00.bl'):
        os.remove(':00.bl')
        
    if removeData:
        if os.path.isfile(polarfile):
            os.remove(polarfile)

        if not saveDatfile:
            if os.path.isfile(defaultDatfile):
                os.remove(defaultDatfile)

    res = {}
    res['cd'] = cd
    res['cl'] = cl
    res['alpha'] = alpha
    res['cm'] = cm
    res['xtr_top'] = xtr_top
    res['xtr_bot'] = xtr_bot
    res['Re'] = Reval
    res['M'] = Mval
    return res

# XFOIL('cl', 
#   'NACA2412',
#   val = 0.2, 
#   minval = 0.0,
#   maxval = 1.0,
#   steps = 10,
#   Re = 1e7,
#   # M = 0.0,
#   flapLocation = None,
#   flapDeflection = 0.0,
#   polarfile='polars.txt',
#   defaultDatfile = 'airfoil.dat',
#   saveDatfile = False,
#   removeData = True,
#   max_iter=100)