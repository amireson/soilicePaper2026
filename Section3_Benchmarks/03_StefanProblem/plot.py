import numpy as np
import matplotlib.pyplot as pl

tL0,zL0=np.loadtxt('Lunardini_0.csv',delimiter=',',skiprows=1,unpack=True)
tL0a,zL0a=np.loadtxt('Lunardini_0a.csv',delimiter=',',skiprows=1,unpack=True)

tL10,zL10=np.loadtxt('Lunardini_10.csv',delimiter=',',skiprows=1,unpack=True)
tL10a,zL10a=np.loadtxt('Lunardini_10a.csv',delimiter=',',skiprows=1,unpack=True)

tL100,zL100=np.loadtxt('Lunardini_100.csv',delimiter=',',skiprows=1,unpack=True)
tL100a,zL100a=np.loadtxt('Lunardini_100a.csv',delimiter=',',skiprows=1,unpack=True)

tNum1,zNum1=np.loadtxt('StefanNum1.csv',delimiter=',',skiprows=1,unpack=True)
tNum2,zNum2=np.loadtxt('StefanNum2.csv',delimiter=',',skiprows=1,unpack=True)
tNum3,zNum3=np.loadtxt('StefanNum3.csv',delimiter=',',skiprows=1,unpack=True)

t,zS1,zN1=np.loadtxt('StefanAna1.csv',delimiter=',',skiprows=1,unpack=True)
t,zS2,zN2=np.loadtxt('StefanAna2.csv',delimiter=',',skiprows=1,unpack=True)
t,zS3,zN3=np.loadtxt('StefanAna3.csv',delimiter=',',skiprows=1,unpack=True)

pl.figure(figsize=(8,4))

pl.subplot(1,2,1)

pl.plot([],[],'k.',markerfacecolor='w',lw=0.75,label='Stefan')
pl.plot([],[],'--k',label='Neumann')

pl.plot(tNum1,zNum1,'-',lw=2,alpha=1.,label=r'$n=0.5$, $T_S=1$, $T_0=-0.01$')
pl.plot(t,zN1,'--k')
pl.plot(t,zS1,'k.',markerfacecolor='w',lw=0.75)


pl.plot(tNum2,zNum2,'-',lw=2,alpha=1.,label=r'$n=0.25$, $T_S=1$, $T_0=-0.01$')
pl.plot(t,zN2,'--k')
pl.plot(t,zS2,'k.',markerfacecolor='w',lw=0.75)

pl.plot(tNum3,zNum3,'-',lw=2,alpha=1.,label=r'$n=0.5$, $T_S=5$, $T_0=-5$')
# pl.plot(t,zS3,'k.',markerfacecolor='w',lw=0.75)
pl.plot(t,zN3,'--k')

# pl.plot(tNum4,zNum4,'-b',alpha=1.,label=r'$n=0.5$, $T_S=5$, $T_0=-5$')
# # pl.plot(t,zS3,'k.',markerfacecolor='w',lw=0.75)
# pl.plot(t,zN4,'--k')

# pl.plot(data.iloc[:,1],'k--',label='Neumann benchmark (Kurylyk et al., 2014)')

# pl.plot(t,m*np.sqrt(t*86400),'-g')

# pl.plot(data.iloc[:,3],'k--')

# pl.plot(tL2,zL2,'-r',alpha=1,label='soilice, v=100m/yr')
# # pl.plot(t,X100,'ko',markersize=4,markerfacecolor='w')
# pl.plot(t,zL2a,'-k.',lw=0.75,label='Lundardini benchmark\n(Kurylyk et al., 2014)')
pl.legend(loc=3,ncols=1,frameon=False)
pl.grid()
pl.ylim(0.34,0)
pl.ylabel('Thaw depth below ground (m)')
pl.xlabel('Time (d)')
pl.title('Solutions with conduction only')


pl.subplot(1,2,2)
pl.plot([],[],'k.',markerfacecolor='w',lw=0.75,label='Stefan')
pl.plot([],[],'--k',label='Lundardini')

pl.plot([],[],'-',lw=2,alpha=1.,label=r'$T_S=1$, $T_0=-0.01$, $q=0$m/yr')
#pl.plot(tL0a,zL0a,'--k')
#pl.plot(t,zS1,'k.',markerfacecolor='w',lw=0.75)

pl.plot(tL10,zL10,'-',lw=2,alpha=1.,label=r'$T_S=1$, $T_0=-0.01$, $q=10$m/yr')
pl.plot(tL10a,zL10a,'--k')

pl.plot(tL100,zL100,'-',lw=2,alpha=1,label=r'$T_S=1$, $T_0=-0.01$, $q=100$m/yr')
pl.plot(tL100a,zL100a,'--k')

pl.plot(tL0,zL0,'-',color='#1f77b4',lw=2,alpha=1.)
pl.plot(tL0a,zL0a,'--k')
pl.plot(t,zS1,'k.',markerfacecolor='w',lw=0.75)



pl.legend(loc=3,frameon=False)
pl.grid()
pl.ylim(0.34,0)
#pl.ylabel('Thaw depth below ground (m)')
pl.xlabel('Time (d)')
pl.title('Solutions including advection')
pl.savefig('thawBenchmarks.png',dpi=300)
