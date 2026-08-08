import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from soilice import model

# from soilice import writeDefaultPars
# writeDefaultPars()

x=model()

x.opts={
 'gravity': 0.0,
 'infiltration': 0.0,
 'cryoK': 0.0,
 'cryoGradient': 0.0,
 'withadv': 0.0,
 'conductionTop': 0.0,
 'conductionBot': 0.0,
 'simulateFlow': 1.0,
 'simulateTransport': 0.0,
 'freeDrainage': 0.0}

# Parameters:
x.readPars()

# Time grid:
x.tGrid(0,1000/24/60.,0.1/24./60.)

# Spatial grid:
x.zGrid(np.arange(0,1.,0.001))

# Boundary conditions:
Se0=0.01
SeI=0.99
psiT=-(SeI**(-1/x.pars['m'])-1)**(1/x.pars['n'])/x.pars['alpha']
psiB=-(Se0**(-1/x.pars['m'])-1)**(1/x.pars['n'])/x.pars['alpha']
print(psiT,psiB)
x.setBCs(psiT=psiT)

# Initial conditions:
x.setICs(T0=1, psi0=psiB)

# Run model
out=x.run()

#to=[5.,10.,20.,50.,100.]
to=[10.,100.,1000.]
ti=[np.where(np.abs(out.t*24*60-i)<1e-4)[0][0] for i in to]
theta=out.thetaT[ti,:]

# Save output
fname='theta_soiliceMathias.csv'
f=open(fname,'w')
f.write('z, psi\n')
for i,j,k,l in zip (out.z,theta[0,:],theta[1,:],theta[2,:]): f.write('%.4f, %.4f, %.4f, %.4f\n'%(i,j,k,l))
f.close()

## Load Mathias solution:
#th,x1,x2,x3,x4,x5=np.loadtxt('%s.csv'%'loam',delimiter=',',unpack=True)
#
#pl.figure(figsize=(5,4))
#
#pl.plot(x1,th,'-k')
#pl.plot(x2,th,'-k')
#pl.plot(x3,th,'-k')
#pl.plot(x4,th,'-k')
#pl.plot(x5,th,'-k')
#pl.plot(out.z*100,theta.T,'.')
#pl.plot(np.nan,np.nan,'-k',label='analytical solution')
#pl.plot(np.nan,np.nan,'.k',label='soilice numerical solution')
#pl.ylim(0.2,0.6)
#pl.xlim(1.,100)
#pl.xscale('log')
#pl.xlabel('Distance (m)',fontsize=13)
#pl.ylabel('Water content (-)',fontsize=13)
#pl.legend(loc=1)
#pl.grid()
#pl.savefig('MathiasProblem.png',dpi=300)
