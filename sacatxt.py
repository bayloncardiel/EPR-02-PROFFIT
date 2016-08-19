import proffitHDF
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy as sp
import scipy.stats
import matplotlib.dates as mdates
import datetime as dt

ff=sys.argv[1]			### ARCHIVO CON FECHAS A UTILIZAR

###########################################################################################################################
### GENERA LISTAS DE HORAS, FECHAS, COLUMNAS, SZAS Y RMS DE UNA SERIE DE FECHAS PARA UN TEST Y UNA ESPECIE DETERMINADA, SE PUEDE ESPECIFICAR HORAS Y RMS

def lee(test,fecha,especie):
	fname = []
	epochTime = []
	columnas = []
	szas = []
	RMSVent = []
	rRMSVent = []
	wShift = []
	TOTerr = []
	#LOSerr = []
	#NOISEerr = []
	for fec in fecha:
		for f_hdf in os.listdir(test+'/'+fec+'/'):
			if os.path.splitext(f_hdf)[1] == '.hdf5':
				try:
					prf = proffitHDF.proffitHDF(test+'/'+fec+'/'+f_hdf)
				except:
					print 'El archivo '+f_hdf+' del dia '+fec+' en '+test+' no se pudo abrir'
					continue
				
				if prf.totCol[prf.especies == especie] > 0:
					VMR_prf = prf.VMR[prf.especies == especie]
					fname.append(f_hdf)
					epochTime.append(prf.epochTime)
					columnas.append(prf.totCol[prf.especies == especie])
					szas.append(prf.SZA)
					RMSVent.append(prf.meanRMS)
					rRMSVent.append(prf.relativeRMS)
					wShift.append(prf.meanWShift)
					TOTerr.append(prf.pcolTOT_err[0])
					#LOSerr.append(np.mean(100.0*prf.LOS_STerr/VMR_prf))
					#NOISEerr.append(np.mean(100.0*prf.NOISE_STerr/VMR_prf))
					#plt.plot(prf.Atot[prf.especies == especie],prf.z)
				else:
					print 'El archivo '+f_hdf+' del dia '+fec+' en '+test+' tiene columna negativa'
	fname = np.array(fname)
	epochTime = np.array(epochTime)
	columnas = np.array(columnas)
	szas = np.array(szas)
	RMSVent = np.array(RMSVent)
	rRMSVent = np.array(rRMSVent)
	wShift = np.array(wShift)
	TOTerr = np.array(TOTerr)
	#LOSerr = np.array(LOSerr)
	#NOISEerr = np.array(NOISEerr)
	
	fname = fname[epochTime.argsort()]
	columnas = columnas[epochTime.argsort()]
	szas = szas[epochTime.argsort()]
	RMSVent = RMSVent[epochTime.argsort()]
	rRMSVent = rRMSVent[epochTime.argsort()]
	wShift = wShift[epochTime.argsort()]
	TOTerr = TOTerr[epochTime.argsort()]
	#LOSerr = LOSerr[epochTime.argsort()]
	#NOISEerr = NOISEerr[epochTime.argsort()]
	epochTime = epochTime[epochTime.argsort()]

	return (fname,epochTime, columnas, szas, RMSVent, rRMSVent, wShift, TOTerr)#LOSerr, NOISEerr)

###########################################################################################################################
###########################################################################################################################
#test1 = "/home/D2_RESULTS/PROFFIT_results/ALTZ/CO2_v1"
#test1 = "/home/D2_RESULTS/PROFFIT_results/ALTZ/O2_v1"
#test1 = "/home/D2_PROFFIT/Jorge_Files/CO2_v1"
#test1 = "/home/D2_PROFFIT/Jorge_Files/O2_v1"
#test1 = "/home/D2_HR125/EM27/RESULTS/O2_EM27_v1"
test1 = "/home/D2_RESULTS/PROFFIT_results/CCA/O2_v1"
gas = test1.split('/')[-1].split('_')[0]
print 'TARGET GAS :',gas
###########################################################################################################################
###########################################################################################################################
print 'PATH:',test1
if os.path.exists(ff):		### COMPRUEBA QUE EL ARCHIVO DE FECHAS EXISTA
	filefechas =open(ff)
	lineas = filefechas.readlines()
	fec = [x.strip('\n') for x in lineas]
	print 'FECHAS', fec
else:	
	print 'NO EXISTE '+ff+' EN '+test1

#fnGAS, etGAS, colGAS, szaGAS, RMSGAS, rRMSGAS, wsGAS, LOSGAS, NOISEGAS = lee(test1,fec,gas)	### LLAMA A LEE()
fnGAS, etGAS, colGAS, szaGAS, RMSGAS, rRMSGAS, wsGAS, TOTGAS = lee(test1,fec,gas)	### LLAMA A LEE()

string = ''
for i,ele in enumerate(fnGAS):
	#string = string + str(ele)+" "+str(etGAS[i])+" %.4f" % szaGAS[i]+" %.4e" % colGAS[i]+" %.4e" % RMSGAS[i]+" %.4e" % rRMSGAS[i]+" %.4e" % wsGAS[i]+" %.6f" % LOSGAS[i]+" %.6f" % NOISEGAS[i]+"\n"
	string = string + str(ele)+" "+str(etGAS[i])+" %.4f" % szaGAS[i]+" %.4e" % colGAS[i]+" %.4e" % RMSGAS[i]+" %.4e" % rRMSGAS[i]+" %.4e" % wsGAS[i]+" %.4e" % TOTGAS[i]+"\n"

fileout = open(min(fec)[0:4]+'-'+max(fec)[0:4]+'_'+gas+'_CCA.txt','w')	### ARCHIVO TXT CON UNA MATRIZ QUE INCLUYE EL NOMBRE DEL ARCHIVO, EL EPOCH TIME, EL SZA, LA COLUMNA, EL RMS Y EL RMS RELATIVO (RMS/SIGNAL) DE CADA MEDICION CORRIDA EN PROFFIT
fileout.write(string)
fileout.close
print 'ARCHIVO '+min(fec)[0:4]+'-'+max(fec)[0:4]+'_'+gas+'_CCA.txt'+' CREADO'
