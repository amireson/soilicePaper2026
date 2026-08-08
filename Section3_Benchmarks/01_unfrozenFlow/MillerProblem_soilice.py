import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from soilice.src_soil import model

# from soilice import writeDefaultPars
# writeDefaultPars()

x=model()

x.opts={'simulateFlow': 1.0,
 'gravity': 1.0,
 'infiltration': 0.0,
 'cryoK': 0.0,
 'cryoGradient': 0.0,
 'freeDrainage': 0.0,
 'simulateTransport': 0.0,
 'withadv': 0.0,
 'conductionTop': 0.0,
 'conductionBot': 0.0}

# Parameters:
x.readPars()

# Time grid:
x.tGrid(0,1.,0.001)

# Spatial grid:
x.zGrid(np.arange(0,4,0.00625))

# Boundary conditions:
psiT=0.1
x.setBCs(psiT=psiT)

# Initial conditions:
x.setICs(T0=1, psi0=[-4,0])

# Run model
out=x.run()

print(out.t[[100,500,1000]])

# Save output
fname='psi_soiliceMiller.csv'
f=open(fname,'w')
f.write('z, psi\n')
for i,j,k,l in zip (out.z,out.psie[100,:],out.psie[500,:],out.psie[1000,:]): f.write('%.4f, %.4f, %.4f, %.4f\n'%(i,j,k,l))
f.close()

