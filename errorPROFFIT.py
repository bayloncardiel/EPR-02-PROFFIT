import proffitHDF
import numpy as np

test_O2 = "/home/D2_RESULTS/PROFFIT_results/ALTZ/O2_v1"
test_CO2 = "/home/D2_RESULTS/PROFFIT_results/ALTZ/CO2_v1"
fec = "151228"
f_hdf = "151228.20_115025SN.bin.hdf5"
prf = proffitHDF.proffitHDF(test_CO2+'/'+fec+'/'+f_hdf)

print np.mean(prf.LOS_STerr)/4.0
print np.mean(prf.NOISE_STerr)/4.0
