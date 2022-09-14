import math as math
import numpy as np

def cosspace(startPoint = -1, endPoint = 1, N = 100):
    theta = np.linspace(0,np.pi,N)
    x = (1-np.cos(theta))/2.
    dpts = endPoint - startPoint
    vals = startPoint + x*dpts
    return vals