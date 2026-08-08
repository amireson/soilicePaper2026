import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
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


#########################
# First results figure
#########################
def nicePlot(run,fn):

    times=[0.1,1,5,10]
    ti=[int(i/run.dt) for i in times]
    # ti=[0, 300]
    lowerlim=0.6

    fig = pl.figure(figsize=(5, 7), constrained_layout=True)

    gs = fig.add_gridspec(
        nrows=3,
        ncols=2,
        height_ratios=[2, 1, 1]
    )

    # Two subplots across the upper half
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Full-width subplots in the lower quarters
    ax3 = fig.add_subplot(gs[1, :])
    ax4 = fig.add_subplot(gs[2, :])

    # Example labels
    ax1.set_title("(a)")
    ax2.set_title("(b)")
    ax3.set_title("(c)")
    ax4.set_title("(d)")

    ########################################################################

    ax1.plot(run.T[ti,:].T,run.z,'r')
    ax1.set_ylim(lowerlim,0)
    # ax1.set_xlim(-0.5,0)
    ax1.grid()
    ax1.set_xlabel(r'Temperature ($^\circ$C)')
    ax1.set_ylabel('Depth (m)')
    ax1.text(-0.9,0.1,r'$t=0.1$ d')
    ax1.text(-0.7,0.2,r'$t=1.0$ d')
    ax1.text(-0.5,0.3,r'$t=5.0$ d')
    ax1.text(-0.3,0.4,r'$t=10.$ d')



    ax2.plot(run.thetaL[ti,:].T-run.pars['thetaR'],run.z,'b')
    ax2.plot(run.thetaI[ti,:].T,run.z,'c')
    ax2.plot(run.thetaT[ti,:].T,run.z,'g')
    ax2.plot([],[],'b',label=r'Liquid, $\theta_l-\theta_r$')
    ax2.plot([],[],'c',label=r'Ice, $\theta_i$')
    ax2.plot([],[],'g',label=r'Total, $\theta_t$')

    ax2.legend(loc=3)
    ax2.set_ylim(lowerlim,0)
    ax2.set_xlim(0,0.55)
    ax2.grid()
    ax2.set_xlabel('Water content (-)')
    ax2.set_yticklabels([])


    t,qT,qB,dm,jT,jB,du,jTA,jBA,duL=getBalance(run)
    # ax3.plot(t[1:],np.cumsum(run.qI[1:]*run.dt)*1000,color='purple',label='Cum. irrigation')
    # ax3.plot(t[1:],np.cumsum(run.qI[1:]*run.dt*1000)-qT[1:],color='blue',label='Cum. runoff')
    # ax3.plot(t,qB,color='red',label='Cum. drainage')
    # ax3.plot(t[1:],qT[1:]-qB[1:],'.',color='orange',label='Net water balance flux')
    # ax3.plot(t[1:],dm[1:],'-k',label='Cumulative change in mass')
    # ax3.grid(); ax3.legend(ncols=2)
    # ax3.set_ylabel(r'Mass (kg m$^{-2}$)')

    ax3.plot(t[1:],run.qI[1:]*1000,'r',label='Irrigation')
    ax3.plot(t[1:],np.diff(qT)/run.dt,'--b',label='Infiltration')
    ax3.plot(t[1:],run.qI[1:]*1000-np.diff(qT)/run.dt,color='purple',label='Runoff')
    ax3.plot(t[1:],np.diff(qB)/run.dt,'--',color='green',label='Drainage')
    ax3.plot(t[1:],np.diff(dm)/run.dt,'-',color='orange',label='Change in mass')
    ax3.plot(t[1:],np.diff(qT-qB)/run.dt,'--k',label='Net flux')
    ax3.grid(); ax3.legend(ncols=3,fontsize=9,frameon=False)
    ax3.set_ylabel('Mass flux \n (kg m$^{-2}$ d$^{-1}$)')

    ax4.plot(t[1:],np.diff(jT-jTA)/1e6/run.dt,'.-',color='green',label='Top conduction')
    # ax4.plot(t[1:],np.diff(jTA)/1e6/run.dt,'--',color='blue',label='Top advection')
    ax4.plot(t[1:],np.diff(jB-jBA)/1e6/run.dt,'-',color='purple',label='Bottom conduction')
    # ax4.plot(t[1:],np.diff(jBA)/1e6/run.dt,'--',color='green',label='Bottom advection')


    ax4.plot(t[1:],np.diff(du)/run.dt/1e6,'-',color='orange',label='Change in internal energy')
    ax4.plot(t[1:],np.diff(jT-jB)/1e6/run.dt,'--k',label='Net flux')
    ax4.grid(); ax4.legend(ncols=2,fontsize=9,frameon=False)
    ax4.set_ylabel('Energy flux \n(MJ m$^{-2}$ d$^{-1}$)')

    pl.savefig(fn,dpi=300)

nicePlot(run_1,'run_1.png')
printBalance(run_1)

#########################
# Second results figure
#########################

def doublePlot(run,old,fn,ncol=3):

    times=[0.1,1,5,10]
    ti=[int(i/run.dt) for i in times]
    ti2=[int(i/old.dt) for i in times]
    lowerlim=0.6
    fig = pl.figure(figsize=(5,7), constrained_layout=True)

    gs = fig.add_gridspec(
        nrows=3,
        ncols=2,
        height_ratios=[2, 1, 1]
    )

    # Two subplots across the upper half
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Full-width subplots in the lower quarters
    ax3 = fig.add_subplot(gs[1, :])
    ax4 = fig.add_subplot(gs[2, :])

    # Example labels
    ax1.set_title("(a)")
    ax2.set_title("(b)")
    ax3.set_title("(c)")
    ax4.set_title("(d)")

    ########################################################################

    ax1.plot(run.T[ti,:].T,run.z,'r')
    ax1.plot(old.T[ti2,:].T,old.z,'r',alpha=0.25)
    ax1.set_ylim(lowerlim,0)
    # ax1.set_xlim(-0.5,0)
    ax1.grid()
    ax1.set_xlabel(r'Temperature ($^\circ$C)')
    ax1.set_ylabel('Depth (m)')
    ax1.text(-0.9,0.1,r'$t=0.1$ d')
    ax1.text(-0.7,0.2,r'$t=1.0$ d')
    ax1.text(-0.5,0.3,r'$t=5.0$ d')
    ax1.text(-0.3,0.4,r'$t=10.$ d')



    ax2.plot(run.thetaL[ti,:].T-run.pars['thetaR'],run.z,'b')
    ax2.plot(old.thetaL[ti2,:].T-old.pars['thetaR'],old.z,'b',alpha=0.25)
    ax2.plot(run.thetaI[ti,:].T,run.z,'c')
    ax2.plot(old.thetaI[ti2,:].T,old.z,'c',alpha=0.25)
    ax2.plot(run.thetaT[ti,:].T,run.z,'g')
    ax2.plot(old.thetaT[ti2,:].T,old.z,'g',alpha=0.25)
    ax2.plot([],[],'b',label=r'Liquid, $\theta_l-\theta_r$')
    ax2.plot([],[],'c',label=r'Ice, $\theta_i$')
    ax2.plot([],[],'g',label=r'Total, $\theta_t$')

    ax2.legend(loc=3)
    ax2.set_ylim(lowerlim,0)
    ax2.set_xlim(0,0.55)
    ax2.grid()
    ax2.set_xlabel('Water content (-)')
    ax2.set_yticklabels([])


    t,qT,qB,dm,jT,jB,du,jTA,jBA,duL=getBalance(run)
    # ax3.plot(t,np.cumsum(run.qI*run.dt)*1000,color='purple',label='Cum. irrigation')
    # ax3.plot(t,np.cumsum(run.qI*run.dt*1000)-qT,color='blue',label='Cum. runoff')
    # ax3.plot(t,qB,color='blue',label='Cum. drainage')
    # ax3.plot(t,qT-qB,'.',color='orange',label='Net water balance')
    # ax3.plot(t,dm,'-k',label='Cum. change in mass')
    # ax3.grid(); ax3.legend(ncol=2)
    # ax3.set_ylabel(r'Mass flux\n (kg m$^{-2}$)')

    ax3.plot(t[1:],run.qI[1:]*1000,'r',label='Irrigation')
    ax3.plot(t[1:],np.diff(qT)/run.dt,'--b',label='Infiltration')
    ax3.plot(t[1:],run.qI[1:]*1000-np.diff(qT)/run.dt,color='purple',label='Runoff')
    ax3.plot(t[1:],np.diff(qB)/run.dt,'--',color='green',label='Drainage')

    ax3.plot(t[1:],np.diff(dm)/run.dt,'-',color='orange',label='Change in mass')
    ax3.plot(t[1:],np.diff(qT-qB)/run.dt,'--k',label='Net flux')
    ax3.grid(); ax3.legend(ncols=ncol,fontsize=9,frameon=False)
    ax3.set_ylabel('Mass flux \n (kg m$^{-2}$ d$^{-1}$)')


    ax4.plot(t[1:],np.diff(jT-jTA)/1e6/run.dt,'.-',color='green',label='Top conduction')
    # ax4.plot(t[1:],np.diff(jTA)/1e6/run.dt,'--',color='blue',label='Top advection')
    ax4.plot(t[1:],np.diff(jB-jBA)/1e6/run.dt,'-',color='purple',label='Bottom conduction')
    # ax4.plot(t[1:],np.diff(jBA)/1e6/run.dt,'--',color='green',label='Bottom advection')


    ax4.plot(t[1:],np.diff(du)/run.dt/1e6,'-',color='orange',label='Change in internal energy')
    ax4.plot(t[1:],np.diff(jT-jB)/1e6/run.dt,'--k',label='Net flux')
    ax4.grid(); ax4.legend(ncols=2,fontsize=9,frameon=False)
    ax4.set_ylabel('Energy flux \n(MJ m$^{-2}$ d$^{-1}$)')

    pl.savefig(fn,dpi=300)

doublePlot(run_2,run_1,'run_2.png',ncol=2)
printBalance(run_2)


#########################
# Third results figure
#########################

doublePlot(run_3,run_1,'run_3.png')
printBalance(run_3)

#########################
# Results 4-7 figure
#########################
def balanceBar(run,x,xlab):
    t,qT,qB,dm,jT,jB,du,jTA,jBA,duL=getBalance(run)
    width=0.9
    JTA=jTA[-1]/1e3
    JTC=(jT[-1]-jTA[-1])/1e3
    JBA=(jBA[-1])/1e3
    JBC=(jB[-1]-jBA[-1])/1e3
    DU=du[-1]/1e3
    DUL=duL[-1]/1e3
    DUS=DU-DUL

    pl.bar(x,JTC,width=width,color='g',label='Top conduction in')
    pl.bar(x,JTA,bottom=JTC,width=width,color='brown',label='Top advection in')
    pl.bar(x,-JBA,bottom=JTA+JTC,width=width,color='r',label='Bottom advection in')

    pl.bar(x+1,JBC,width=width,color='b',label='Bottom conduction out')
    pl.bar(x+1,DUS,bottom=JBC,width=width,color='c',hatch='///',label='Increase in Sensible heat')

    if DUL<0:
        pl.bar(x,-100,bottom=-100,width=width,color='orange',hatch='\\\\\\',label=r'$\Delta$ Latent heat - thawing')
        pl.bar(x,-DUL,bottom=JTA+JTC-JBA,width=width,color='skyblue',hatch='\\\\\\',label=r'$\Delta$ Latent heat - freezing')

    else:
        pl.bar(x+1,DUL,bottom=DUS+JBC,width=width,color='orange',hatch='\\\\\\',label=r'$\Delta$ Latent heat - thawing')
        pl.bar(x+1,-100,bottom=-100,width=width,color='skyblue',hatch='\\\\\\',label=r'$\Delta$ Latent heat - freezing')

    # pl.text(1,3960,'+ve',ha='center')
    # pl.text(2,3960,'-ve',ha='center')
    pl.xlabel(xlab) 
    pl.gca().set_xticks([])
    pl.ylim(0,4000)




pl.figure(figsize=(10,9))
ax1=pl.subplot(2,4,1)
pl.plot(run_4.T[-1,:],run_4.z,'k',label='Run 4')
pl.plot(run_5.T[-1,:],run_4.z,color='orange',label='Run 5')
pl.plot(run_6.T[-1,:],run_4.z,'--r',label='Run 6')
pl.plot(run_7.T[-1,:],run_4.z,'--g',label='Run 7')
pl.plot(run_4.T[0,:],run_4.z,'--b',label='Initial condition')
pl.ylim(0.6,0)
pl.grid()
pl.xlabel(r'$T$ ($^\circ$C)')
pl.ylabel('Depth (m)')
ax1.set_title("(a)")

ax2=pl.subplot(2,4,2)
pl.plot(run_4.thetaI[-1,:],run_4.z,'k',label='Run 4')
pl.plot(run_5.thetaI[-1,:],run_4.z,color='orange',label='Run 5')
pl.plot(run_6.thetaI[-1,:],run_4.z,'--r',label='Run 6')
pl.plot(run_7.thetaI[-1,:],run_4.z,'--g',label='Run 7')
pl.plot(run_4.thetaI[0,:],run_4.z,'--b',label='Initial condition')
pl.ylim(0.6,0)
pl.grid()
pl.gca().set_yticklabels([])
pl.xlabel(r'$\theta_i$ (-)')
ax2.set_title("(b)")

ax3=pl.subplot(2,4,3)
pl.plot(run_4.thetaL[-1,:],run_4.z,'k',label='Run 4: \nbaseline')
pl.plot(run_5.thetaL[-1,:],run_4.z,color='orange',label='Run 5: \nno infiltration')
pl.plot(run_6.thetaL[-1,:],run_4.z,'--r',label='Run 6: \nno advection')
pl.plot(run_7.thetaL[-1,:],run_4.z,'--g',label='Run 7: \nno surface conduction')
pl.plot(run_4.thetaL[0,:],run_4.z,'--b',label='Initial condition')
pl.ylim(0.6,0)
pl.grid()
pl.legend(ncols=5,loc=4,bbox_to_anchor=(2.15, -0.35))
pl.gca().set_yticklabels([])
pl.xlabel(r'$\theta_l$ (-)')
ax3.set_title("(c)")

ax4=pl.subplot(2,4,4)
pl.plot(run_4.thetaT[-1,:],run_4.z,'k',label='Run 4: baseline')
pl.plot(run_5.thetaT[-1,:],run_4.z,color='orange',label='Run 5: no infiltration')
pl.plot(run_6.thetaT[-1,:],run_4.z,'--r',label='Run 6')
pl.plot(run_7.thetaT[-1,:],run_4.z,'--g',label='Run 7')
pl.plot(run_4.thetaT[0,:],run_4.z,'--b',label='Initial condition')
pl.ylim(0.6,0)
pl.grid()
pl.gca().set_yticklabels([])
pl.xlabel(r'$\theta_t$ (-)')
# pl.subplots_adjust(wspace=0.05)
ax4.set_title("(d)")

# pl.figure(figsize=(10,4))
ax5=pl.subplot(2,4,5)
balanceBar(run_4,1,'Run 4: Baseline')
pl.ylabel('Energy exchanged (kJ m$^{-2}$)')
ax5.set_title("(e)")

ax6=pl.subplot(2,4,6)
balanceBar(run_5,1,'Run 5: No infiltration')
pl.gca().set_yticklabels([])
ax6.set_title("(f)")

ax7=pl.subplot(2,4,7)
balanceBar(run_6,1,'Run 6: No advection')
pl.gca().set_yticklabels([])
pl.legend(ncols=4,loc=4,bbox_to_anchor=(2.25, -0.33),fontsize=10)
ax7.set_title("(g)")

ax8=pl.subplot(2,4,8)
balanceBar(run_7,1,'Run 7: No surface conduction')
pl.gca().set_yticklabels([])
ax8.set_title("(h)")

# pl.text(-500,1.5,'(Refreezing)',fontsize=10)

# pl.legend(framealpha=1,loc=2)
pl.subplots_adjust(wspace=0.05,hspace=0.47)

pl.savefig('run4-7.png',dpi=300)


print('\nRUN4:'); printBalance(run_4)
print('\nRUN5:'); printBalance(run_5)
print('\nRUN6:'); printBalance(run_6)
print('\nRUN7:'); printBalance(run_7)
