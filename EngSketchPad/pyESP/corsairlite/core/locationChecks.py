def inIPython():
    """
    Checks to see if user is operating in an iPython notebook for the 
    purposes of displaing output.
    """
    try:
        cfg = get_ipython().config 
        if 'IPKernelApp' in cfg.keys():
            if 'connection_file' in cfg['IPKernelApp'].keys():
                return False
            else:
                return True
        else:
            return True
    except NameError:
        return False

def inNotebook():
    """
    Checks to see if user is operating in an iPython notebook for the 
    purposes of displaing output.
    """
    try:
        cfg = get_ipython().config 
        if 'IPKernelApp' in cfg.keys():
            if 'connection_file' in cfg['IPKernelApp'].keys():
                return True
            else:
                return False
        else:
            return False
    except NameError:
        return False