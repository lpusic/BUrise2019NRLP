#from brian2 import reinit, clear, PoissonGroup, Hz, SpikeMonitor, run, seed, ms, second
import brian2 as brian
from brian2 import magic_network
import brian2.core.magic
import numpy as np
import os, sys

nruns = int(sys.argv[1])

for nrun in range(1, nruns+1):
    brian.seed(nrun)
    print('RUN: ' + str(nrun))
    #brian2.core.magic.reinit(states = True)
    #brian.clear(erase = True, all = True)
    rate = int(sys.argv[2])
    foldername = 'rate'+str(rate)+'/run_'+str(nrun)
    os.system('mkdir -p -v '+foldername)
    
    N  = 1000
    time_input = 23000 * brian.ms
    P  = brian.PoissonGroup(N,rate * brian.Hz)
    S = brian.SpikeMonitor(P)
    
    P.rates = rate * brian.Hz
    brian.run(time_input, report='text', report_period = 10 * brian.second)
    
    fname = 'noise_' 
    #print(type(S.spike_trains()[0][0]/brian.ms))   
    for s in range(len(S.spike_trains().keys())):
        spiketimes = [round(1000*x/brian.ms,1)+50 for x in list(S.spike_trains()[s])]
        np.savetxt(foldername+'/'+fname+str(s)+'.txt',spiketimes,fmt='%10.1f',newline='\n')
    
    

