import numpy as np
import matplotlib.pyplot as pl
from soilice import model
from soilice.src_constitutiveFunctions import (
    thetaFun, CFun, KFun,
    thermalKfun, GCEFun, SFCslope, CBFun
    )
from soilice import writeDefaultPars
from numba import types
from numba.typed import Dict

def MakeDictFloat():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,)
    return d

x=model()
x.readPars()
pars=MakeDictFloat()
const=MakeDictFloat()
for i in x.pars: pars[i]=x.pars[i]
for i in x.const: const[i]=x.const[i]

I=10.
E=0.

n=200
log=False

fig=pl.figure(figsize=(7,5))

if log: 
    psie=np.logspace(-2,2,n)
    xx=-psie
else:
    psie=np.linspace(-10,1,n)
    xx=psie

Tarr=np.array([0,-.005,-.01,-0.05,-1])

pl.subplot(2,2,1)
for T in Tarr:
    psif=np.full(n,GCEFun(np.array([T]),pars,const))
    pars['E']=0.
    pars['impedance']=0.
    pl.plot(xx,KFun(np.minimum(psif,psie),np.minimum(psif,psie),pars,const),label=fr'$T = $ {T}')
pl.grid()
#pl.ylabel('$K$ (m d$^{-1}$)]\n as a function of liquid VWC')
pl.ylabel('$K$ (m d$^{-1}$)]')
# pl.xlabel(r'$T$ (deg C)')
if log: pl.gca().set_xscale('log')
if log: pl.xlim(1,1e-4)
pl.gca().set_yscale('log')
# yl=pl.gca().get_ylim()
pl.ylim(1e-6,1)
#pl.gca().set_xticklabels([])
pl.legend()
pl.xlabel(r'$\psi_e$ (m)')
pl.title(r'$K = f(\psi_f)$ (Equation 50)')

pl.subplot(2,2,2)
for T in Tarr:
    psif=np.full(n,GCEFun(np.array([T]),pars,const))
    pars['E']=E
    pars['impedance']=I
    pl.plot(xx,KFun(np.full(n,psie),np.minimum(psif,psie),pars,const),label=str(T))
pl.grid()
pl.gca().set_yscale('log')
#pl.ylabel('$K$ (m d$^{-1}$)]\#n with impedance model')

pl.xlabel(r'$\psi_e$ (m)')
if log: pl.gca().set_xscale('log')
if log: pl.xlim(1,1e-4)
# pl.ylim(0.005,0.022)
pl.ylim(1e-6,1)

if log: 
    T=-np.logspace(-4,0,n)
    psif=GCEFun(T,pars,const)
else:
    T=np.linspace(-1,0.2,n)
    psif=GCEFun(T,pars,const)
    T=-T
psi_arr=np.array([0,-.5,-1,-2,-5])
pl.gca().set_yticklabels([])

pl.title(r'$K = f(\psi_e,\theta_i)$ (Equation 51)')

pl.subplot(2,2,3)
for psie in psi_arr:
    pars['impedance']=0.
    pars['E']=0.
    pl.plot(-T,KFun(np.minimum(psif,psie),np.minimum(psif,psie),pars,const),label=fr'$\psi_e = $ {psie}')
pl.grid()
# pl.xlabel(r'$T$ (deg C)')
if log: pl.gca().set_xscale('log')
if log: pl.xlim(1,1e-4)
pl.gca().set_yscale('log')
# yl=pl.gca().get_ylim()
pl.ylim(1e-6,1)
#pl.gca().set_xticklabels([])
# pl.gca().set_yticklabels([])
pl.legend()
pl.xlabel(r'$T$ (deg C)')
pl.ylabel('$K$ (m d$^{-1}$)]')

pl.subplot(2,2,4)
for psie in psi_arr:
    pars['impedance']=I
    pars['E']=E
    pl.plot(-T,KFun(np.full(n,psie),np.minimum(psif,psie),pars,const),label=str(psie))
pl.grid()
pl.gca().set_yscale('log')
pl.xlabel(r'$T$ (deg C)')
if log: pl.gca().set_xscale('log')
if log: pl.xlim(1,1e-4)
# pl.ylim(0.005,0.022)
pl.ylim(1e-6,1)
pl.gca().set_yticklabels([])


pl.subplots_adjust(hspace=0.3,wspace=0.07,top=0.93)
pl.savefig('KFuns.png',dpi=300)
