import numpy as np
import matplotlib.pyplot as pl
from scipy.special import erfc

from soilice import model
from soilice.src_soil import MakeDictFloat
from soilice.src_constitutiveFunctions import thermalKfun, CBFun

def OgataBanks(t,x,c0,D,v):
    A=erfc((x-v*t)/2/np.sqrt(D*t))+np.exp(v*x/D)*erfc((x+v*t)/2/np.sqrt(D*t))
    c=c0*0.5*A
    return c
    
s=model()

s.opts={
 'gravity': 1.0,
 'infiltration': 1.0,
 'cryoK': 0.0,
 'cryoGradient': 0.0,
 'withadv': 1.0,
 'conductionTop': 1.0,
 'conductionBot': 0.0,
 'simulateFlow': 1.0,
 'simulateTransport': 1.0,
 'freeDrainage': 1.0}

s.tGrid(0,20,2)
s.zGrid(np.arange(0,4.01,0.01))

s.readPars()

psi0=-0.41655092025858204 # Obtained for SS flow
s.setICs(T0=0.,psi0=psi0)
s.setBCs(qI=0.1,TInf=1.,TTop=1.)

o=s.run()

print('Water content is ',o.thetaL[-1,-1],o.thetaT[0,0])
np.array([o.psie[-1,-1]])

parsD = MakeDictFloat() 
for k in s.pars: parsD[k] = s.pars[k]

constD = MakeDictFloat()
for k in s.const: constD[k] = s.const[k]

kappa=thermalKfun(np.array([o.psie[-1,-1]]),np.array([o.psie[-1,-1]]),np.array([1]),parsD,constD)
CB=CBFun(np.array([o.psie[-1,-1]]),np.array([o.psie[-1,-1]]),parsD,constD)
D=kappa/CB
v=s.qI[0]*s.const['rho_liq']*s.const['cp_liq']/CB

print('Kappa: ',kappa/86400)
print('CB:    ',CB)
print('D:     ',D)
print('v:     ',v)

# i=np.arange(20,s.nt,20)
pl.figure(figsize=(5,4))
pl.plot(o.T[1:,:].T,s.z,'-',color='#1f77b4',lw=2)
for ti in s.t[1:]:
    Ta=OgataBanks(ti,s.z,s.TTop[0],D,v)
    pl.plot(Ta,s.z,'k--',markersize=4,markerfacecolor='w')
pl.ylim(2,0)
pl.grid()
pl.plot(0,3,'-',color='#1f77b4',lw=2,label='soilice')
pl.plot(0,3,'k--',label='Ogata Banks')
pl.legend()
pl.ylabel('Depth (m)')
pl.xlabel(r'Temperature ($^o$C)')
pl.savefig('ogatabanks.png',dpi=300)
pl.close(pl.gcf())

#pl.figure(figsize=(7,7))
#o.plotBalance()
#pl.show()

