import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from soilice import model
from soilice import save

# Set parameters:
TI=0.
TG=0.
T0=-1.
qI=3/1000.
psi0=-1.

x=model()
x.opts={'gravity': 1.0,
 'infiltration': 1.0,
 'cryoK': 0.0,
 'cryoGradient': 0.0,
 'withadv': 1.0,
 'conductionTop': 1.0,
 'conductionBot': 1.0,
 'simulateFlow': 1.0,
 'simulateTransport': 1.0,
 'freeDrainage': 1.0}

dz=0.002
zMax=1.0
x.zGrid(np.arange(0,zMax+dz,dz))

x.tGrid(0,10,0.01)

x.readPars()
x.pars['impedance']=10.

x.setBCs(qI=qI,TInf=TI,TBot=T0,TTop = TG)
x.setICs(psi0=psi0,T0=T0)

out=x.run()

save('out.dill',out)