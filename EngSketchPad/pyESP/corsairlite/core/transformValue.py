import numpy as np

def transform(rawval, transform, inversionFlag = False):
    if transform is not None:
        if not inversionFlag:
            if transform in ['log','ln']:
                tval = np.log(rawval)
            if transform == 'log10':
                tval = np.log10(rawval)
            if transform == 'exp':
                tval = np.exp(rawval)
            if transform in ['10^', '10**']:
                tval = 10.0**rawval
            if transform[0] == '+':
                tval = rawval + float(transform[1:])
            if transform[0] == '-':
                tval = rawval - float(transform[1:])
            if transform[0] == '*' and transform[0:2] != '**':
                tval = rawval * float(transform[1:])
            if transform[0] == '/':
                tval = rawval / float(transform[1:])
            if transform[0] == '^' or transform[0:2] == '**':
                tval = rawval ** float(transform[2:])

        else:
            if transform in ['log','ln']:
                tval = np.exp(rawval)
            if transform == 'log10':
                tval = 10.0**rawval
            if transform == 'exp':
                tval = np.log(rawval)
            if transform in ['10^', '10**']:
                tval = np.log10(rawval)
            if transform[0] == '+':
                tval = rawval - float(transform[1:])
            if transform[0] == '-':
                tval = rawval + float(transform[1:])
            if transform[0] == '*' and transform[0:2] != '**':
                tval = rawval / float(transform[1:])
            if transform[0] == '/':
                tval = rawval * float(transform[1:])
            if transform[0] == '^' or transform[0:2] == '**':
                tval = rawval ** (1.0/float(transform[2:]))

        return tval
    else:
        return rawval