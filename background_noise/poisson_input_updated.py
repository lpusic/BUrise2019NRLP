import brian2
import numpy as np
import os, sys

nruns = int(sys.argv[1])

for nrun in range(1, nruns+1):
    brian2.seed(nrun)
    print ('RUN: ' + str(nrun))
    brian2.start_scope()
    rate = int(sys.argv[2])
    foldername = 'rate'+str(rate)+'/run_'+str(nrun)
    os.system('mkdir -p -v '+foldername)
    
    N  = 1000
    time_input = 23000 * brian2.ms
    P  = brian2.PoissonGroup(N, rate * brian2.Hz)
    S  = brian2.SpikeMonitor(P)

    brian2.run(time_input, report='stdout', 
               report_period = 2 * brian2.second)
    
    fname = 'noise_'
    spiketimes_all = S.spike_trains()
    for s in range(len(spiketimes_all)):
        spiketimes = [round(1000*x / brian2.second,1)+50 for x in list(spiketimes_all[s])]
        np.savetxt(foldername+'/'+fname+str(s)+'.txt',spiketimes,fmt='%.1f',newline='\n')
    
    


