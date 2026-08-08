import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

pl.figure(figsize=(8,4))
pl.subplot(1,2,1)
fname='psi_CeliaMiller.csv'
z,psi1,psi2,psi3=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)

fname='psi_soiliceMiller.csv'
zs,psi1s,psi2s,psi3s=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)

ti=[100,500,1000]
pl.plot(psi1s,zs,lw=2,label='soilice; t=0.1 d')
pl.plot(psi2s,zs,lw=2,label='soilice; t=0.5 d')
pl.plot(psi3s,zs,lw=2,label='soilice; t=1.0 d')
pl.plot([],[],'--k',label='Celia solution')
pl.plot(psi1,z,'--k')
pl.plot(psi2,z,'--k')
pl.plot(psi3,z,'--k')
pl.ylim(4,0)

#
pl.grid()
pl.xlabel('Pressure head (m)')
pl.ylabel('Depth (m)')
pl.legend(loc=3)
#pl.subplots_adjust(top=0.96)
#pl.savefig('Miller.png',dpi=300)
pl.title('Ponded Infiltration Problem\nMiller et al. (1998)')

pl.subplot(1,2,2)
# Load Mathias solution:
th,x1,x2,x3=np.loadtxt('theta_Mathias.csv',delimiter=',',unpack=True)

# Load soilice solution:
fname='theta_soiliceMathias.csv'
x,th1,th2,th3=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)

pl.plot(x,th1,'-',lw=2,label='soilice\nt=10min')
pl.plot(x,th2,'-',lw=2,label='soilice\nt=100min')
pl.plot(x,th3,'-',lw=2,label='soilice\nt=1000min')
pl.plot(x1/100,th,'--k')
pl.plot(x2/100,th,'--k')
pl.plot(x3/100,th,'--k')
#
#pl.plot(x1/100,th,'-k')
#pl.plot(x2/100,th,'-k')
#pl.plot(x3/100,th,'-k')
#pl.plot(x4/100,th,'-k')
#pl.plot(x5/100,th,'-k')
#pl.plot(x,th1,'.')
#pl.plot(x,th2,'.')
#pl.plot(x,th3,'.')
#pl.plot(x,th4,'.')
#pl.plot(x,th5,'.')
#pl.plot(out.z*100,theta.T,'.')
pl.plot(np.nan,np.nan,'--k',label='analytical\nsolution')
pl.ylim(0.2,0.55)
pl.xlim(.001,1)
pl.xscale('log')
pl.xlabel('Distance (m)')
pl.ylabel('Water content (-)')
pl.legend(loc=3)
pl.grid()
pl.title('Horizontal infiltration\nMathias & Sander (2021)')

pl.subplots_adjust(wspace=0.35)
pl.savefig('unfrozenFlow.png',dpi=300)
