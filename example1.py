"""
This code reproduces example 1 from the paper
"Optimal monotonicity-preserving perturbations of a given Runge-Kutta method".
In the paper results are given for just a few methods; here we include
several others.
"""
from nodepy import rk
import numpy as np

#beta = 0.5#3-2*sqrt(2)

def f(t,u):
    if u<0:
        raise Exception('Negative value encountered')
    elif u-1>1.e-15:
        raise Exception('Value greater than 1 encountered')
    else:
        return np.sign(np.sin(t))*u*(1.-u)
        #return (9.*np.sign(np.sin(t))+7)*u*(1.-u)*(beta-u)
        #return (9.*np.sign(np.sin(t))-7)*u*(1.-u)*(beta-u)
        #return 2.*np.sign(np.sin(t))*u*(1.-u)*(beta-u)
        #return 2.*(np.sin(t))*u*(1.-u)*(beta-u)

from nodepy.ivp import IVP

ly_problem = IVP(f=f, u0=0.49, T=100.)


methods = []
methods.append(rk.loadRKM('FE'))
methods.append(rk.loadRKM('Mid22'))
methods.append(rk.loadRKM('MTE22'))
methods.append(rk.loadRKM('SSP22'))
methods.append(rk.loadRKM('SSP22star'))
methods.append(rk.loadRKM('Heun33'))
methods.append(rk.loadRKM('SSP33'))
methods.append(rk.loadRKM('RK44'))
methods.append(rk.loadRKM('Merson43'))
methods.append(rk.loadRKM('SSP104'))
methods.append(rk.loadRKM('Fehlberg45'))
methods.append(rk.loadRKM('DP5'))
methods.append(rk.loadRKM('BS5'))
methods.append(rk.loadRKM('SSP75'))
methods.append(rk.loadRKM('SSP85'))
methods.append(rk.loadRKM('SSP95'))
methods.append(rk.loadRKM('CMR6'))
methods.append(rk.loadRKM('PD8'))

for method in methods:
    hmax = 7.1
    hmin = 0.001
    eps = 0.001
    while hmax-hmin>eps:
        h = 0.5*(hmax+hmin)
        failed = False
        for u0 in [1.e-8,1-1.e-8]:
            try:
                ly_problem.u0 = u0
                t,y = method(ly_problem,dt=h)
            except:
                failed = True
            if failed:
                hmax = h
                break
        if not failed:
            hmin = h

    #print 'done'
    print(method.shortname)
    print("h_obs = "+str(hmin))
    print('')
