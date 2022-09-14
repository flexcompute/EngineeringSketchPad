import struct
import copy
import numpy as np

def readSensx(filename = 'sensx.afl'):
    f = open(filename,'rb')
    allData = f.read()
    f.close()

    checksumArray = []

    buffer = allData

    sizeof_int = 4
    sizeof_char = 1
    sizeof_double = 8

    pointer = 0

    def movePointer(buffer, pointer, moveBytes, typeKey, notrunc=False):
        dta = struct.unpack(typeKey,buffer[pointer:pointer+moveBytes])
        pointer += moveBytes
        if notrunc:
            return dta, pointer
        else:
            return dta[0], pointer
        
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    dta,           pointer = movePointer(buffer,pointer,32,'32s')           # MSES
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    dta,           pointer = movePointer(buffer,pointer,32,'32s')           # Airfoil
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    kalpha,        pointer = movePointer(buffer,pointer,sizeof_int,'i')
    kmach,         pointer = movePointer(buffer,pointer,sizeof_int,'i')
    kreyn,         pointer = movePointer(buffer,pointer,sizeof_int,'i')
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    ldepma,        pointer = movePointer(buffer,pointer,sizeof_int,'i')
    ldepre,        pointer = movePointer(buffer,pointer,sizeof_int,'i')
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    alpha,         pointer = movePointer(buffer,pointer,sizeof_double,'d')
    mach,          pointer = movePointer(buffer,pointer,sizeof_double,'d')
    reyn,          pointer = movePointer(buffer,pointer,sizeof_double,'d')
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    dnrms,         pointer = movePointer(buffer,pointer,sizeof_double,'d')
    drrms,         pointer = movePointer(buffer,pointer,sizeof_double,'d')
    dvrms,         pointer = movePointer(buffer,pointer,sizeof_double,'d')
    dnmax,         pointer = movePointer(buffer,pointer,sizeof_double,'d')
    drmax,         pointer = movePointer(buffer,pointer,sizeof_double,'d')
    dvmax,         pointer = movePointer(buffer,pointer,sizeof_double,'d')
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    ii,            pointer = movePointer(buffer,pointer,sizeof_int,'i')
    nbl,           pointer = movePointer(buffer,pointer,sizeof_int,'i')
    nmod,          pointer = movePointer(buffer,pointer,sizeof_int,'i')
    npos,          pointer = movePointer(buffer,pointer,sizeof_int,'i')
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    ileb = []
    iteb = []
    for i in range(0,nbl):
        dta,           pointer = movePointer(buffer,pointer,sizeof_int,'i')
        ileb.append(dta)
        dta,           pointer = movePointer(buffer,pointer,sizeof_int,'i')
        iteb.append(dta)
    end, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    xleb   = []
    yleb   = []
    xteb   = []
    yteb   = []
    sblegn = []
    for i in range(0,nbl):
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        xleb.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        yleb.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        xteb.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        yteb.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        sblegn.append(dta)
    end, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
    checksumArray.append(beg == end)
    beg,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
    cl,            pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cm,            pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdw,           pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdv,           pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdf,           pointer = movePointer(buffer,pointer,sizeof_double,'d')
    al_alfa,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cl_alfa,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cm_alfa,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdw_alfa,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdv_alfa,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdf_alfa,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    al_mach,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cl_mach,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cm_mach,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdw_mach,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdv_mach,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdf_mach,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    al_reyn,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cl_reyn,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cm_reyn,       pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdw_reyn,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdv_reyn,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    cdf_reyn,      pointer = movePointer(buffer,pointer,sizeof_double,'d')
    end,           pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
        
    iend = [0]*len(ileb)
    for i in range(0,nbl):
        ioff = ileb[i]-1
        ileb[i] -= ioff
        iteb[i] -= ioff
        iend[i] = ii-ioff

    xbi = []
    ybi = []
    cp = []
    hk = []
    cp_alfa = []
    hk_alfa = []
    cp_mach = []
    hk_mach = []
    cp_reyn = []
    hk_reyn = []

    for i in range(0,nbl):
        j = iend[i] - ileb[i]
        for k in range(0,2):
            xbi.append(np.zeros(j))
            ybi.append(np.zeros(j))
            cp.append(np.zeros(j))
            hk.append(np.zeros(j))
            cp_alfa.append(np.zeros(j))
            hk_alfa.append(np.zeros(j))
            cp_mach.append(np.zeros(j))
            hk_mach.append(np.zeros(j))
            cp_reyn.append(np.zeros(j))
            hk_reyn.append(np.zeros(j))
            beg, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
            for m in range(0,j):
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                xbi[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                ybi[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                cp[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                hk[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                cp_alfa[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                hk_alfa[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                cp_mach[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                hk_mach[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                cp_reyn[2*i+k][m] = dta
                dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                hk_reyn[2*i+k][m] = dta
            end, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
            checksumArray.append(beg == end)
            

    modn    = []
    al_mod  = []
    cl_mod  = []
    cm_mod  = []
    cdw_mod = []
    cdv_mod = []
    cdf_mod = []
    gn = [[None for jjj in range(0,2*(nbl-1)+2)] for i in range(0,nmod)]
    xbi_mod = [[None for jjj in range(0,2*(nbl-1)+2)] for i in range(0,nmod)]
    ybi_mod = [[None for jjj in range(0,2*(nbl-1)+2)] for i in range(0,nmod)]
    cp_mod = [[None for jjj in range(0,2*(nbl-1)+2)] for i in range(0,nmod)]
    hk_mod = [[None for jjj in range(0,2*(nbl-1)+2)] for i in range(0,nmod)]

    # for k in range(0,nmod):   
    for k in range(0,40):   
        beg, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        modn.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        al_mod.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        cl_mod.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        cm_mod.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        cdw_mod.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        cdv_mod.append(dta)
        dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
        cdf_mod.append(dta)
        end, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
        checksumArray.append(beg == end)

        for ib in range(0,nbl):
            for ix in range(0,2):
                j = iteb[ib] - ileb[ib] + 1
                beg, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
                dta, pointer = movePointer(buffer,pointer,sizeof_double*j,'d'*j, notrunc=True)
                gn[k][2*ib+ix] = dta
                end, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
                checksumArray.append(beg == end)
                
                j = iend[ib] - ileb[ib]
                beg, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &begMarker
                xbi_row = []
                ybi_row = []
                cp_row = []
                hk_row = []
                for m in range(0,j):
                    dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                    xbi_row.append(dta)
                    dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                    ybi_row.append(dta)
                    dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                    cp_row.append(dta)
                    dta, pointer = movePointer(buffer,pointer,sizeof_double,'d')
                    hk_row.append(dta)
                end, pointer = movePointer(buffer,pointer,sizeof_int,'i')     # &endMarker
                checksumArray.append(beg == end)
                
                xbi_mod[k][2*ib+ix] = xbi_row
                ybi_mod[k][2*ib+ix] = ybi_row
                cp_mod[k][2*ib+ix]  = cp_row
                hk_mod[k][2*ib+ix]  = hk_row

    if all(checksumArray):
        # print('Markers all match')
        pass
    else:
        raise ValueError('Import markers do not match')

    outDict = {}

    outDict['kalpha']    = kalpha
    outDict['kmach']     = kmach
    outDict['kreyn']     = kreyn
    outDict['ldepma']    = ldepma
    outDict['ldepre']    = ldepre
    outDict['alpha']     = alpha
    outDict['mach']      = mach
    outDict['reyn']      = reyn
    outDict['dnrms']     = dnrms
    outDict['drrms']     = drrms
    outDict['dvrms']     = dvrms
    outDict['dnmax']     = dnmax
    outDict['drmax']     = drmax
    outDict['dvmax']     = dvmax
    outDict['ii']        = ii
    outDict['nbl']       = nbl
    outDict['nmod']      = nmod
    outDict['npos']      = npos
    outDict['ileb']      = ileb
    outDict['iteb']      = iteb
    outDict['xleb']      = xleb
    outDict['yleb']      = yleb
    outDict['xteb']      = xteb
    outDict['yteb']      = yteb
    outDict['sblegn']    = sblegn
    outDict['cl']        = cl
    outDict['cm']        = cm
    outDict['cdw']       = cdw
    outDict['cdv']       = cdv
    outDict['cdf']       = cdf
    outDict['al_alfa']   = al_alfa
    outDict['cl_alfa']   = cl_alfa
    outDict['cm_alfa']   = cm_alfa
    outDict['cdw_alfa']  = cdw_alfa
    outDict['cdv_alfa']  = cdv_alfa
    outDict['cdf_alfa']  = cdf_alfa
    outDict['al_mach']   = al_mach
    outDict['cl_mach']   = cl_mach
    outDict['cm_mach']   = cm_mach
    outDict['cdw_mach']  = cdw_mach
    outDict['cdv_mach']  = cdv_mach
    outDict['cdf_mach']  = cdf_mach
    outDict['al_reyn']   = al_reyn
    outDict['cl_reyn']   = cl_reyn
    outDict['cm_reyn']   = cm_reyn
    outDict['cdw_reyn']  = cdw_reyn
    outDict['cdv_reyn']  = cdv_reyn
    outDict['cdf_reyn']  = cdf_reyn
    outDict['xbi']       = xbi
    outDict['ybi']       = ybi
    outDict['cp']        = cp
    outDict['hk']        = hk
    outDict['cp_alfa']   = cp_alfa
    outDict['hk_alfa']   = hk_alfa
    outDict['cp_mach']   = cp_mach
    outDict['hk_mach']   = hk_mach
    outDict['cp_reyn']   = cp_reyn
    outDict['hk_reyn']   = hk_reyn
    outDict['modn']      = modn
    outDict['al_mod']    = al_mod
    outDict['cl_mod']    = cl_mod
    outDict['cm_mod']    = cm_mod
    outDict['cdw_mod']   = cdw_mod
    outDict['cdv_mod']   = cdv_mod
    outDict['cdf_mod']   = cdf_mod
    outDict['gn']        = gn
    outDict['xbi_mod']   = xbi_mod
    outDict['ybi_mod']   = ybi_mod
    outDict['cp_mod']    = cp_mod
    outDict['hk_mod']    = hk_mod


    return outDict