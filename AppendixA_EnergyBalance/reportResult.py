import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from soilice import loadModel

run_1=loadModel('run_1.dill')
run_0=loadModel('../Section4_Infiltration/run_1.dill')

def getBalance(self):

    t=self.t
    nt,nz=self.thetaL.shape
    ml=self.thetaL*self.dz*self.const['rho_liq']
    mi=self.thetaI*self.dz*self.const['rho_ice']
    ms=np.zeros((nt,nz))+((1-self.pars['thetaS'])*self.dz*self.pars['rho_soil'])
    u=(ml*self.const['cp_liq']+mi*self.const['cp_ice']+ms*self.pars['cp_soil'])*self.T-mi*self.           const['lambda_f']
    ml=np.sum(ml,axis=1)
    mi=np.sum(mi,axis=1)
    u=np.sum(u,axis=1)
    du=u-u[0]
    m=ml+mi
    dm=m-m[0]
    qT=self.qT.cumsum()*self.const['rho_liq']
    qB=self.qB.cumsum()*self.const['rho_liq']
    jT=self.jT.cumsum()
    jB=self.jB.cumsum()
    jTA=qT*self.const['cp_liq']*self.TInf*self.opts['withadv']
    jBA=qB*self.const['cp_liq']*self.T[:,-1]*self.opts['withadv']
    uL=-mi*self.const['lambda_f']
    duL=uL-uL[0]
    return t,qT,qB,dm,jT,jB,du,jTA,jBA,duL

def printBalance(run):
    t,qT,qB,dm,jT,jB,du,jTA,jBA,duL,=getBalance(run)
    print(f'Total mass irrigated:        {run.qI[1:].sum()*1000*run.dt:.2f} kg/m2')
    print(f'Total mass infiltrated:      {qT[-1]:.2f} kg/m2')
    print(f'Total mass of runoff:        {run.qI[1:].sum()*1000*run.dt-qT[-1]:.2f} kg/m2')
    print(f'Total mass of drainage:      {qB[-1]:.2f} kg/m2')
    print(f'Change in mass:              {dm[-1]:.2f} kg/m2')
    print(f'Mass balance error:          {dm[-1]-(qT[-1]-qB[-1]):.4f} kg/m2')
    print(f'MB error as fraction of qT:  {(dm[-1]-(qT[-1]-qB[-1]))/qT[-1]:.4f} -')
    print('')
    print(f'Total top advection:         {jTA[-1]/1e3:.2f} kJ/m2')
    print(f'Total top conduction:        {(jT[-1]-jTA[-1])/1e3:.2f} kJ/m2')
    print(f'Total bottom advection:      {(jBA[-1])/1e3:.2f} kJ/m2')
    print(f'Total bottom conduction:     {(jB[-1]-jBA[-1])/1e3:.2f} kJ/m2')
    print(f'Change in internal energy:   {du[-1]/1e3:.1f} kJ/m2')
    print(f'Change in Latent heat:       {duL[-1]/1e3:.1f} kJ/m2')
    print(f'Change in Sensible heat:     {(du[-1]-duL[-1])/1e3:.1f} kJ/m2')
    print(f'Energy balance error:        {(du[-1]-(jT[-1]-jB[-1]))/1e3:.4f} kJ/m2')
    print(f'EB error as fraction of jT:  {(du[-1]-(jT[-1]-jB[-1]))/jT[-1]:.4f} -')


print('Original model run:')
printBalance(run_0)

print('\nSimplified MB model run:')
printBalance(run_1)

def RMSE(x,y):
    return np.sqrt(np.mean((x-y)**2))

print('')
print(f'RMSE of thetaT is {RMSE(np.squeeze(run_1.thetaT),np.squeeze(run_0.thetaT))}')
print(f'RMSE of thetaL is {RMSE(np.squeeze(run_1.thetaL),np.squeeze(run_0.thetaL))}')
print(f'RMSE of thetaI is {RMSE(np.squeeze(run_1.thetaI),np.squeeze(run_0.thetaI))}')
print(f'RMSE of Temperature is {RMSE(np.squeeze(run_1.T),np.squeeze(run_0.T))}')

