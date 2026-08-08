#######################################################
# 
# Analytical solutions from Kurylyk, B. L., McKenzie, J. M., MacQuarrie, K. T. B., 
# and Voss, C. I. (2014). Analytical solutions for benchmarking cold regions
# subsurface water flow and energy transport models: One-dimensional soil thaw with
# conduction and advection. Advances in Water Resources, 70:172–184.
#
# Used to validate the soilice model in Ireson et al. 2026 (in-prep)
#
# Coded up by A. Ireson, 28th July 2026
#
#######################################################

import numpy as np
from scipy.optimize import fsolve
from scipy.special import erf, erfc

# Needed constituitive functions:
def thermalKfun(psie,psif,T,pars,const):
    # Uses geometric mean of each component
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    thetaA=pars['thetaS']-thetaT
    kappa = (
        (const['kappa_liq'] ** thetaL) *
        (const['kappa_ice'] ** thetaI) *
        (const['kappa_air'] ** thetaA) *
        (pars['kappa_soil'] ** pars['theta_mineral'])*
        (pars['kappa_org'] ** pars['theta_org']))
    return kappa

def thetaFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

def CBFun(thetaL,thetaI,pars,const):
    thetaS=1-pars['thetaS']
    c_pb=(const['cp_ice']*thetaI*const['rho_ice'])+(const['cp_liq']*thetaL*const['rho_liq'])+(pars['cp_soil']*thetaS*pars['rho_soil'])
    return c_pb

def get_th_kappa_cpb(psi,pars,const):
    psi=np.array([psi])
    thetaL=thetaFun(psi,pars)   # Total water content = liquid water content in thawed zone

    # Thawed properties:    
    kappa=thermalKfun(psi,psi,psi,pars,const)
    c_pb=CBFun(thetaL,0,pars,const)

    # Frozen properties:
    kappa_f=thermalKfun(psi,np.array([-1e6]),psi,pars,const)
    c_pb_f=CBFun(0,thetaL,pars,const)
    
    return thetaL,kappa,c_pb,kappa_f,c_pb_f

# For the Stefan and Lunardini solutions - function to evaluate ST
def StefanNumber(Ts,Tf,thetaL,c_pb,const):
    ST=c_pb*(Ts-Tf) /const['rho_liq'] /const['lambda_f'] /thetaL
    return ST

# The Stefan solution for the depth of freezing, X
def StefanSoln(t,Ts,Tf,thetaL,c_pb,kappa,const):
    ST=StefanNumber(Ts,Tf,thetaL,c_pb,const)
    alpha=kappa/c_pb
    X=np.sqrt(2*ST*alpha*t)
    return X

# For the Neumann solution - function to evaluate Y1-Y2
def fN(m,Ts,Tf,Ti,c_u,c_f,kappa_u,kappa_f,thetaL,const):
    Y1=0.5*const['lambda_f']*thetaL*const['rho_liq']*np.sqrt(np.pi)*m
    Y2left=(
        np.sqrt(c_u*kappa_u)*(Ts-Tf)*
        (np.exp(-m**2*c_u/4/kappa_u)/
        erf(m/2/np.sqrt(kappa_u/c_u)))
    )
    Y2right=(
        np.sqrt(c_f*kappa_f)*(Tf-Ti)*
        (np.exp(-m**2*c_f/4/kappa_f)/
         erfc(m/2/np.sqrt(kappa_f/c_f)))
    )
    Y2=(Y2left-Y2right)
    # print(Y1,Y2left,Y2right,Y1-Y2)
    return Y1-Y2

# The Neumann solution for the depth of freezing, X
# Coded up in the form provided in Ireson et al, 2026 in-prep
def Neumman(t,Ts,Tf,Ti,alpha_u,alpha_f,kappa_u,kappa_f,thetaL,const):
    c_u=kappa_u/alpha_u
    c_f=kappa_f/alpha_f
    m0=0.00001
    m=fsolve(fN, x0=m0, args=(Ts,Tf,Ti,c_u,c_f,kappa_u,kappa_f,thetaL,const))[0]
    X=m*np.sqrt(t*86400)
    # print(f'm = {m}')
    return X,m

# For the Lunardini solution - function to solve to equal zero
def f(X, vt, St, t, alpha, Vscale):
    # Lunardini's solution equalling zero
    V=vt*Vscale
    return X + alpha/V*(np.exp(-V*X/alpha) - 1) - V*St*t

# The Lunardini solution for the depth of freezing, X
def LunardiniSoln(t,Ts,Tf,vt,thetaL,c_pb,kappa,const):
    alpha=kappa/c_pb
    ST=StefanNumber(Ts,Tf,thetaL,c_pb,const)
    Vscale=const['rho_liq']*const['cp_liq']/c_pb
    X0 = np.array([np.sqrt(2*ST*alpha*ti) for ti in t])  # Stefan solution as first guess for iteration
    X=np.array([fsolve(f, x0=X0[i], args=(vt, ST, t[i], alpha, Vscale)) for i in range(len(t))])
    return X
