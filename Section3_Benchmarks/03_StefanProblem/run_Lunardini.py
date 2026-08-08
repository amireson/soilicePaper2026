import numpy as np
import matplotlib.pyplot as pl

from analyticalSoln import StefanSoln
from analyticalSoln import LunardiniSoln
from analyticalSoln import get_th_kappa_cpb

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
    # Kurylyk linear SFC model:
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
    o=x.run()
    print('Energy input = ',o.jT.sum())
    
    # Get the zero degree isotherm depth from the soilice numerical model output
    fig, ax = pl.subplots()
    cs = ax.contour(o.t, o.z, o.T.T, [0])
    segments = cs.allsegs[0]
    pl.close(fig)
    all_points = np.vstack(segments)
    t_vals = all_points[:, 0]
    z_vals = all_points[:, 1]
    return t_vals,z_vals,x.pars,x.const

# Run  soilice with q=~0
tL0,zL0,parsL0,constL0=runModel(dt=0.01,dx=0.001,Ts=1,T0=-0.01,psi=0,porosity=0.50,q=0.00001)
# Save output
fname='Lunardini_0.csv'
f=open(fname,'w')
f.write('t, X\n')
for i,j in zip (tL0,zL0): f.write('%.4f, %.4f\n'%(i,j))
f.close()


# Run  soilice with q=10/365 mm/d
tL10,zL10,parsL10,constL10=runModel(dt=0.01,dx=0.001,Ts=1,T0=-0.01,psi=0,porosity=0.50,q=10/365)
# Save output
fname='Lunardini_10.csv'
f=open(fname,'w')
f.write('t, X\n')
for i,j in zip (tL10,zL10): f.write('%.4f, %.4f\n'%(i,j))
f.close()


# Run  soilice with q=100/365 mm/d
tL100,zL100,parsL100,constL100=runModel(dt=0.01,dx=0.001,Ts=1,T0=-0.01,psi=0,porosity=0.50,q=100/365)
fname='Lunardini_100.csv'
f=open(fname,'w')
f.write('t, X\n')
for i,j in zip (tL100,zL100): f.write('%.4f, %.4f\n'%(i,j))
f.close()

# Run analytical solutions
t=np.logspace(-2,np.log10(20),40)

thetaL=parsL10['thetaS']
kappa=parsL10['kappau']
c_pb=parsL10['cu']
zL0a=np.squeeze(LunardiniSoln(t,1,-0.001,parsL0['q'],thetaL,c_pb,kappa,constL0))

fname='Lunardini_0a.csv'
f=open(fname,'w')
f.write('t, X\n')
for i,j in zip (t,zL0a): f.write('%.4f, %.4f\n'%(i,j))
f.close()


thetaL=parsL10['thetaS']
kappa=parsL10['kappau']
c_pb=parsL10['cu']
zL10a=np.squeeze(LunardiniSoln(t,1,-0.001,parsL10['q'],thetaL,c_pb,kappa,constL10))

fname='Lunardini_10a.csv'
f=open(fname,'w')
f.write('t, X\n')
for i,j in zip (t,zL10a): f.write('%.4f, %.4f\n'%(i,j))
f.close()

thetaL=parsL100['thetaS']
kappa=parsL100['kappau']
c_pb=parsL100['cu']
zL100a=np.squeeze(LunardiniSoln(t,1,-0.001,parsL100['q'],thetaL,c_pb,kappa,constL100))

fname='Lunardini_100a.csv'
f=open(fname,'w')
f.write('t, X\n')
for i,j in zip (t,zL100a): f.write('%.4f, %.4f\n'%(i,j))
f.close()


