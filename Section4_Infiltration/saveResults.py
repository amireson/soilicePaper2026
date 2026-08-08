import numpy as np
import pandas as pd
from soilice import loadModel

run_1=loadModel('run_1.dill')
run_2=loadModel('run_2.dill')
run_3=loadModel('run_3.dill')
run_4=loadModel('run_4.dill')
run_5=loadModel('run_5.dill')
run_6=loadModel('run_6.dill')
run_7=loadModel('run_7.dill')

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

# Save output from each run in a csv file
def saveOutput(run,fn):

    # Variables with depth, for final timestep
    f=open(f'{fn}_z.csv','w')
    f.write('Depth (m), T (deg C), thetaI (-), thetaL (-), thetaT (-)\n')
    for i in range(run.nz):
        f.write(f'    {run.z[i]:.3f},   {run.T[-1,i]:.4f},')
        f.write(f'     {run.thetaI[-1,i]:.4f},')
        f.write(f'     {run.thetaL[-1,i]:.4f},')
        f.write(f'     {run.thetaT[-1,i]:.4f}\n')
    f.close()

    # Variables with time
    t,qT,qB,dm,jT,jB,du,jTA,jBA,duL=getBalance(run)
    qT=np.diff(qT)/run.dt
    qB=np.diff(qB)/run.dt
    dS=np.diff(dm)/run.dt
    jT=np.diff(jT)/run.dt/1000
    jTA=np.diff(jTA)/run.dt/1000
    jTC=jT-jTA
    jB=np.diff(jB)/run.dt/1000
    jBA=np.diff(jBA)/run.dt/1000
    jBC=jB-jBA
    du=np.diff(du)/run.dt/1000
    duL=np.diff(duL)/run.dt/1000
    duS=du-duL
    runoff=np.maximum(0,run.qI[1:]*1000-qT)

    f=open(f'{fn}_MB.csv','w')
    f.write('Time (d), qI (mm/d), qT (mm/d), qR (mm/d), qB (mm/d), dS (mm/d)\n')
    for i in range(1,run.nt):
        f.write(f'   {t[i]:5.2f},')
        f.write(f'     {run.qI[i]*1000:5.3f},')
        f.write(f'     {qT[i-1]:5.3f},')
        f.write(f'     {runoff[i-1]:5.3f},')
        f.write(f'     {qB[i-1]:5.3f},')
        f.write(f'     {dS[i-1]:5.3f}')
        f.write('\n')
    f.close()

    f=open(f'{fn}_EB.csv','w')
    f.write('Time (d),      jTA (*),      jTC (*),      jBA (*),      jBC (*),      duL (*),      duS (*) * all units are kJ/m2/d\n')
    for i in range(1,run.nt):
        f.write(f'   {t[i]:5.2f},')
        f.write(f'     {jTA[i-1]:8.2f},')
        f.write(f'     {jTC[i-1]:8.2f},')
        f.write(f'     {jBA[i-1]:8.2f},')
        f.write(f'     {jBC[i-1]:8.2f},')
        f.write(f'     {duL[i-1]:8.2f},')
        f.write(f'     {duS[i-1]:8.2f},')
        f.write('\n')
    f.close()

saveOutput(run_1,'run_1')
saveOutput(run_2,'run_2')
saveOutput(run_3,'run_3')
saveOutput(run_4,'run_4')
saveOutput(run_5,'run_5')
saveOutput(run_6,'run_6')
saveOutput(run_7,'run_7')


