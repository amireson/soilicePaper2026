# Benchmarking analytical thawing solutions

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from analyticalSoln import StefanSoln
from analyticalSoln import get_th_kappa_cpb
from analyticalSoln import Neumman

from soilice import model

# Function to run  soilice for analytical benchmarks:
def runModel(dt,dx,Ts,T0,psi,porosity,q=0,Tm=-0.0005,rtol=1e-7):
    x=model(rtol=rtol)
    x.opts={'gravity': 0.0,
     'infiltration': 1.0,
     'cryoK': 0.0,
     'cryoGradient': 0.0,
     'withadv': 1.0,
     'conductionTop': 1.0,
     'conductionBot': 0.0,
     'simulateFlow': 0.0,
     'simulateTransport': 1.0,
     'freeDrainage': 0.0}
    x.zGrid(np.arange(0,2.000001,dx))
    # x.zGrid(np.logspace(-3,np.log10(1.5),200)-0.001)
    x.tGrid(0,20.,dt)
    x.setBCs(TTop=Ts)
    x.setICs(T0=T0,psi0=psi)

    # Kurylyk parameters
    x.readPars()
    pars=x.pars
    pars['lambda_f']=334000.
    pars['rho_ice']=1000.
    pars['rho_liq']=1000.
    pars['thetaS']=porosity
    pars['thetaR']=0.
    pars['Tm']=Tm
    pars['q']=q

    if porosity==0.5:
        pars['kappau']=85867.44095410302
        pars["kappaf"]=170194.561024182
        pars["cu"]=3195000.0
        pars["cf"]=2155000.0
    if porosity==0.25:
        pars["kappau"]=5789.660231777036
        pars["kappaf"]=8151.014321722309
        pars["cu"]=2702500.0
        pars["cf"]=2182500.0

    x.pars=pars
    const=x.const
    const['rho_ice']=pars['rho_ice']
    x.const=const
    o=x.run()
    print('Energy input = ',o.jT.sum())

    # Get the zero degree isotherm depth from the soilice numerical model output
    fig, ax = pl.subplots()
    cs = ax.contour(o.t, o.z, o.T.T, [0])
    segments = cs.allsegs[0]
    pl.close(fig)
    all_points = np.vstack(segments)
    tNum = all_points[:, 0]
    zNum = all_points[:, 1]

    # Get Stefan solution:
    t=np.logspace(-2,np.log10(x.t[-1]),40)
    zS=StefanSoln(t,Ts,T0,porosity,pars['cu'],pars['kappau'],x.const)

    #Get Neumann solution:
    alpha=pars['kappau']/pars['cu']/86400
    alpha_f=pars['kappaf']/pars['cf']/86400
    Tf=0.
    zN,m=Neumman(t,Ts,Tf,T0,alpha,alpha_f,pars['kappau']/86400,pars['kappaf']/86400,porosity,x.const)

    print(zN)
    return tNum,zNum,t,zS,zN

psi=0.

porosity=0.5
Ts=1.
T0=-0.01
tNum1,zNum1,t,zS1,zN1=runModel(dt=0.01,dx=0.001,Ts=Ts,T0=T0,psi=psi,porosity=porosity,Tm=-0.01)

porosity=0.25
Ts=1
T0=-0.01
tNum2,zNum2,t,zS2,zN2=runModel(dt=0.01,dx=0.001,Ts=Ts,T0=T0,psi=psi,porosity=porosity,Tm=-0.01)

porosity=0.5
Ts=5.
T0=-5.
tNum3,zNum3,t,zS3,zN3=runModel(dt=0.01,dx=0.001,Ts=Ts,T0=T0,psi=psi,porosity=porosity,Tm=-0.01)

# Save outputs:
f=open('StefanNum1.csv','w')
f.write('t, X\n')
for i,j in zip (tNum1,zNum1): f.write('%.4f, %.4f\n'%(i,j))
f.close()
f=open('StefanNum2.csv','w')
f.write('t, X\n')
for i,j in zip (tNum2,zNum2): f.write('%.4f, %.4f\n'%(i,j))
f.close()
f=open('StefanNum3.csv','w')
f.write('t, X\n')
for i,j in zip (tNum3,zNum3): f.write('%.4f, %.4f\n'%(i,j))
f.close()

f=open('StefanAna1.csv','w')
f.write('t, S, N\n')
for i,j,k in zip (t,zS1,zN1): f.write('%.4f, %.4f, %.4f\n'%(i,j,k))
f.close()
f=open('StefanAna2.csv','w')
f.write('t, S, N\n')
for i,j,k in zip (t,zS2,zN2): f.write('%.4f, %.4f, %.4f\n'%(i,j,k))
f.close()
f=open('StefanAna3.csv','w')
f.write('t, S, N\n')
for i,j,k in zip (t,zS3,zN3): f.write('%.4f, %.4f, %.4f\n'%(i,j,k))
f.close()
