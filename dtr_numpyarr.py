import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import time
import glob
import proffitHDF

######################################################################
####################	DECIMAL YEAR CALCULATION
def decimal_year(dattim):
	
	def sinceEpoch(date): # returns seconds since epoch
		return time.mktime(date.timetuple())
	s = sinceEpoch

	dec_yr = []
	for ele in dattim:
    		year = ele.year
		startOfThisYear = dt.datetime(year=year, month=1, day=1)
		startOfNextYear = dt.datetime(year=year+1, month=1, day=1)

		yearElapsed = s(ele) - s(startOfThisYear)
		yearDuration = s(startOfNextYear) - s(startOfThisYear)
		fraction = yearElapsed/yearDuration
		dec_yr.append(ele.year + fraction)
	return np.array(dec_yr)

######################################################################
####################	FOURIER SERIES
def fourier_series(diff_d,meas_val,species):
	K = np.ones(shape=(len(diff_d),6))
	for ii,ele in enumerate(K[:,0]):
		K[ii,1] = diff_d[ii]
		K[ii,2] = np.cos(2.0*np.pi*diff_d[ii])
		K[ii,3] = np.sin(2.0*np.pi*diff_d[ii])
		K[ii,4] = np.cos(4.0*np.pi*diff_d[ii])
		K[ii,5] = np.sin(4.0*np.pi*diff_d[ii])

	KT = np.transpose(K)
	mat = np.dot(np.linalg.inv(np.dot(KT,K)),KT)	
	coef = np.dot(mat,meas_val)
	print 'For %s (%i): a0 = %.4e, alpha = %.4e, a1 = %.4e, a2 = %.4e, b1 = %.4e, b2 = %.4e' % (species,len(meas_val),coef[0],coef[1],coef[2],coef[4],coef[3],coef[5])
	#fourier_fun = np.zeros(len(XCO2_unq_d))
	fourier_fun = np.zeros(len(diff_d))
	trend = np.array([coef[0] + coef[1]*ele for ele in diff_d])
	for ii,ele in enumerate(K[:,0]):
		for jj,cof in enumerate(coef):
			fourier_fun[ii] = fourier_fun[ii]+K[ii,jj]*cof
	return fourier_fun,trend,coef

######################################################################
path_dir = '/home/D2_RESULTS/PROFFIT_results/ALTZ/'
method = 'CO2_v2'
especie = 'CO2'
fecha = '150426'
lista_files = glob.glob(path_dir+method+'/'+fecha+'/*.hdf5')

fname = []
epochTime = []
columnas = []
szas = []
for fec in lista_files:
	try:
		prf = proffitHDF.proffitHDF(fec)
	except:
		print 'El archivo %s en %s no se pudo abrir' % (fec.split('/')[-1],path_dir+method+'/'+fecha)
		continue
				
	if prf.totCol[prf.especies == especie] > 0:
		VMR_prf = prf.VMR[prf.especies == especie]
		fname.append(fec)
		epochTime.append(prf.epochTime)
		columnas.append(prf.totCol[prf.especies == especie][0])
		szas.append(prf.SZA)
	else:
		print 'El archivo %s en %s tiene columna negativa' % (fec.split('/')[-1],path_dir+method+'/'+fecha)

fname = np.array(fname)
epochTime = np.array(epochTime)
columnas = np.array(columnas)
szas = np.array(szas)
idx_sort = epochTime.argsort()
fname = fname[idx_sort]
epochTime = epochTime[idx_sort]
columnas = columnas[idx_sort]

datetime_CO2 = [dt.datetime.fromtimestamp(ele) for ele in epochTime]
diff_d_CO2 = np.array(decimal_year(datetime_CO2))	
fourier_fun_CO2,trend_CO2,coef_CO2 = fourier_series(diff_d_CO2,columnas,especie)
dtr_ff_CO2 = np.subtract(columnas,fourier_fun_CO2)

dset_dtr_XCO2 = np.empty(len(epochTime),dtype=[('filenamehdf',np.str_,100),('tepoch',int),('totcol',float)])
dset_dtr_XCO2['filenamehdf'] = fname
dset_dtr_XCO2['tepoch'] = epochTime
dset_dtr_XCO2['totcol'] = dtr_ff_CO2

plt.figure(1)
plt.plot(epochTime,fourier_fun_CO2,linestyle='--')
plt.plot(epochTime,columnas,linestyle='None',marker='+')
plt.figure(2)
plt.plot(epochTime,dtr_ff_CO2,linestyle='None',marker='+')
plt.show()
np.save('/home/D2_RESULTS/PROFFIT_results/TIMESERIES_analysis/npdata/dtrCO2_v3_QA.npy',dset_dtr_XCO2)







