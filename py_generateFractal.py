from neurodsp import sim as ndsp
import numpy as np

# cycle through samples
vals = np.concatenate((np.arange(0.1,1,0.1),np.arange(1,10,1),np.arange(10,110,10)),axis=0)
for samp in np.arange(0,np.size(vals)):
    
    # generate 100 time-series of [samp] duration
    ts = np.zeros((int(vals[samp]*100),1000))
    for i in np.arange(0,1000):
        ts[:,i] = ndsp.aperiodic.sim_powerlaw(vals[samp],100,frange=[5,25])
        
    # write to csv
    np.savetxt('E:/bjg335/projects/reinstatement_fidelity/data/supp_fractal/sig_'+str(samp+1)+'.csv',ts,delimiter=',')