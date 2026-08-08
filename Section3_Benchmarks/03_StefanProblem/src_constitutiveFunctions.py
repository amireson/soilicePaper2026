import numpy as np
from numba import njit

# Editted functions:
@njit(inline='always')
def thetaFun(psi,pars):
    psim=pars['Tm']*334560/(9.81*273.15)
    Se=np.zeros(len(psi))
    Se[psi>psim]=(psi[psi>psim]-psim)/(0-psim)
    Se[psi>0]=1.
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

@njit(inline='always')
def CFun(psi,pars):
    psim=pars['Tm']*334560/(9.81*273.15)
    C=np.zeros(len(psi))
    C[psi>psim]=(pars['thetaS']-pars['thetaR'])/(0-psim)
    C[psi>=0]=0.
    return C

# Thermal conductivity function
@njit(inline='always')
def thermalKfun(psie,psif,T,pars,const):
    psim=pars['Tm']*334560/(9.81*273.15)
    Se=np.zeros(len(psie))
    Se[psif>psim]=(psif[psif>psim]-psim)/(0-psim)
    Se[psif>0]=1.
    kappa=Se*(pars['kappau']-pars['kappaf'])+pars['kappaf']
    return kappa

# Bulk heat capacity function
@njit(inline='always')
def CBFun(psie,psif,pars,const):
    psim=pars['Tm']*334560/(9.81*273.15)
    Se=np.zeros(len(psie))
    Se[psif>psim]=(psif[psif>psim]-psim)/(0-psim)
    Se[psif>0]=1.
    CB=Se*(pars['cu']-pars['cf'])+pars['cf']
    return CB

# Original functions: 
@njit(inline='always')
def KFun(psie,psif,pars,const):
    # Impedance model for K after Taylor and Luthin
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)

    Se=(1+(psie*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psie>0.]=1.0
    Ke=pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

    K=Ke*10**(-10*thetaI)
    return K

# Slope function dtheta_L/dT
@njit(inline='always')
def SFCslope(psie,psif,pars,const):
    C=const['lambda_f']/const['g']/const['T0']
    dthdpsi=CFun(psif,pars)
    dthdT=C*dthdpsi
    dthdT[psie==psif]=0.
    return dthdT
    
# GCE
@njit(inline='always')
def GCEFun(T,pars,const):
    psi=T*const['lambda_f']/(const['g']*const['T0'])
    psi[psi>0]=0.
    return psi
