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

#fechas =['140215','151216']			### ARCHIVO CON FECHAS A UTILIZAR
fechas = ['150418','140428']
bms = ['CaF2','KBr']
lsty = ['-','--']

###########################################################################################################################
### GENERA LISTAS DE HORAS, FECHAS, COLUMNAS, SZAS Y RMS DE UNA SERIE DE FECHAS PARA UN TEST Y UNA ESPECIE DETERMINADA, SE PUEDE ESPECIFICAR HORAS Y RMS


def lee(test,fecha,especie):
	fname = []
	epochTime = []
	columnas = []
	szas = []
	BASE_STerr = []
	ILS_STerr = []
	LOS_STerr = []
	SOLAR_STerr = []
	T_STerr = []
	SPEC_STerr = []
	NOISE_STerr = []
	BASE_SYerr = []
	ILS_SYerr = []
	LOS_SYerr = []
	SOLAR_SYerr = []
	T_SYerr = []
	SPEC_SYerr = []
	NOISE_SYerr = []
	TOTAL_STerr = []

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
				BASE_STerr.append(np.mean(100.0*prf.BASE_STerr/VMR_prf))
				ILS_STerr.append(np.mean(100.0*prf.ILS_STerr/VMR_prf))
				LOS_STerr.append(np.mean(100.0*prf.LOS_STerr/VMR_prf))
				SOLAR_STerr.append(np.mean(100.0*prf.SOLAR_STerr/VMR_prf))
				T_STerr.append(np.mean(100.0*prf.T_STerr/VMR_prf))
				SPEC_STerr.append(np.mean(100.0*prf.SPEC_STerr/VMR_prf))
				NOISE_STerr.append(np.mean(100.0*prf.NOISE_STerr/VMR_prf))
				BASE_SYerr.append(np.mean(100.0*prf.BASE_SYerr/VMR_prf))
				ILS_SYerr.append(np.mean(100.0*prf.ILS_SYerr/VMR_prf))
				LOS_SYerr.append(np.mean(100.0*prf.LOS_SYerr/VMR_prf))
				SOLAR_SYerr.append(np.mean(100.0*prf.SOLAR_SYerr/VMR_prf))
				T_SYerr.append(np.mean(100.0*prf.T_SYerr/VMR_prf))
				SPEC_SYerr.append(np.mean(100.0*prf.SPEC_SYerr/VMR_prf))
				NOISE_SYerr.append(np.mean(100.0*prf.NOISE_SYerr/VMR_prf))

				#plt.plot(prf.Atot[prf.especies == especie],prf.z)
			else:
				print 'El archivo '+f_hdf+' del dia '+fec+' en '+test+' tiene columna negativa'
	fname = np.array(fname)
	epochTime = np.array(epochTime)
	columnas = np.array(columnas)
	szas = np.array(szas)
	BASE_STerr = np.array(BASE_STerr)
	ILS_STerr = np.array(ILS_STerr)
	LOS_STerr = np.array(LOS_STerr)
	SOLAR_STerr = np.array(SOLAR_STerr)
	T_STerr = np.array(T_STerr)
	SPEC_STerr = np.array(SPEC_STerr)
	NOISE_STerr = np.array(NOISE_STerr)
	BASE_SYerr = np.array(BASE_SYerr)
	ILS_SYerr = np.array(ILS_SYerr)
	LOS_SYerr = np.array(LOS_SYerr)
	SOLAR_SYerr = np.array(SOLAR_SYerr)
	T_SYerr = np.array(T_SYerr)
	SPEC_SYerr = np.array(SPEC_SYerr)
	NOISE_SYerr = np.array(NOISE_SYerr)

	
	sorted_idx = epochTime.argsort()
	fname = fname[sorted_idx]
	columnas = columnas[sorted_idx]
	szas = szas[sorted_idx]
	BASE_STerr = BASE_STerr[sorted_idx]
	ILS_STerr = ILS_STerr[sorted_idx]
	LOS_STerr = LOS_STerr[sorted_idx]
	SOLAR_STerr = SOLAR_STerr[sorted_idx]
	T_STerr = T_STerr[sorted_idx]
	SPEC_STerr = SPEC_STerr[sorted_idx]
	NOISE_STerr = NOISE_STerr[sorted_idx]
	BASE_SYerr = BASE_SYerr[sorted_idx]
	ILS_SYerr = ILS_SYerr[sorted_idx]
	LOS_SYerr = LOS_SYerr[sorted_idx]
	SOLAR_SYerr = SOLAR_SYerr[sorted_idx]
	T_SYerr = T_SYerr[sorted_idx]
	SPEC_SYerr = SPEC_SYerr[sorted_idx]
	NOISE_SYerr = NOISE_SYerr[sorted_idx]
	epochTime = epochTime[sorted_idx]

	err_matrix = np.column_stack((fname,epochTime,szas,columnas,BASE_STerr,ILS_STerr,LOS_STerr,SOLAR_STerr,T_STerr,SPEC_STerr,NOISE_STerr,	\
	BASE_SYerr,ILS_SYerr,LOS_SYerr,SOLAR_SYerr,T_SYerr,SPEC_SYerr,NOISE_SYerr))

	return err_matrix

###########################################################################################################################
###########################################################################################################################
testCO2 = "/home/D2_RESULTS/PROFFIT_results/ALTZ/CO2_vERR"
testO2 = "/home/D2_RESULTS/PROFFIT_results/ALTZ/O2_vERR"
###########################################################################################################################
###########################################################################################################################
for ii,fec in enumerate(fechas):
	if (os.path.exists(fec+'_CO2_err.txt') & os.path.exists(fec+'_O2_err.txt')):
		dt_txt = np.dtype([('fname',np.str_,30),('epochTime',int),('SZA',float),('totCol',float), \
		('BASE_ST',float),('ILS_ST',float),('LOS_ST',float),('SOLAR_ST',float),('T_ST',float), \
		('SPEC_ST',float),('NOISE_ST',float),('BASE_SY',float),('ILS_SY',float),('LOS_SY',float),('SOLAR_SY',float), \
		('T_SY',float),('SPEC_SY',float),('NOISE_SY',float)])
		err_matrixCO2 = np.genfromtxt(fec+'_CO2_err.txt',delimiter=' ',dtype=dt_txt)
		err_matrixO2 = np.genfromtxt(fec+'_O2_err.txt',delimiter=' ',dtype=dt_txt)
		try:
			tm_CO2 = np.array([dt.datetime.fromtimestamp(ele).replace(2001,01,01) for ele in err_matrixCO2['epochTime']])
			tm_O2 = np.array([dt.datetime.fromtimestamp(ele).replace(2001,01,01) for ele in err_matrixO2['epochTime']])

			lw_n = 2.8

			plt.figure(1)
			plt.suptitle('CO2 STATISTICAL ERRORS')
			plt.yscale('log')
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['BASE_ST'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(err_matrixCO2['SZA'],err_matrixCO2['ILS_ST'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['LOS_ST'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['SOLAR_ST'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['T_ST'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			#plt.plot(err_matrixCO2['SZA'],err_matrixCO2['SPEC_ST'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['NOISE_ST'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')		
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('CO2_ST.eps')
		
			plt.figure(2)
			plt.suptitle('CO2 SYSTEMATIC ERRORS')
			plt.yscale('log')
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['BASE_SY'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(err_matrixCO2['SZA'],err_matrixCO2['ILS_SY'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['LOS_SY'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['SOLAR_SY'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['T_SY'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			plt.plot(err_matrixCO2['SZA'],err_matrixCO2['SPEC_SY'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			#plt.plot(err_matrixCO2['SZA'],err_matrixCO2['NOISE_SY'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('CO2_SY.eps')

			plt.figure(3)
			plt.suptitle('O2 STATISTICAL ERRORS')
			plt.yscale('log')
			plt.plot(err_matrixO2['SZA'],err_matrixO2['BASE_ST'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(err_matrixO2['SZA'],err_matrixO2['ILS_ST'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['LOS_ST'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['SOLAR_ST'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['T_ST'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			#plt.plot(err_matrixO2['SZA'],err_matrixO2['SPEC_ST'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['NOISE_ST'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('O2_ST.eps')

			plt.figure(4)
			plt.suptitle('O2 SYSTEMATIC ERRORS')
			plt.yscale('log')
			plt.plot(err_matrixO2['SZA'],err_matrixO2['BASE_SY'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(err_matrixO2['SZA'],err_matrixO2['ILS_SY'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['LOS_SY'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['SOLAR_SY'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['T_SY'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			plt.plot(err_matrixO2['SZA'],err_matrixO2['SPEC_SY'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			#plt.plot(err_matrixO2['SZA'],err_matrixO2['NOISE_SY'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('O2_SY.eps')

		except:
			print 'FNAME',err_matrixCO2['fname'],'SZA',err_matrixCO2['SZA']
			#print 'CO2: BASE ST = %.06f, ILS ST = %.06e, LOS ST = %.06f, SOLAR ST = %.06f, T ST = %.06f, SPEC ST = %.06f, NOISE ST = %.06f, BASE SY = %.06f, ILS SY = %.06e, LOS SY = %.06f, SOLAR SY = %.06f, T SY = %.06f, SPEC SY = %.06f, NOISE SY = %.06f' %(err_matrixCO2['BASE_ST'],err_matrixCO2['ILS_ST'],err_matrixCO2['LOS_ST'],err_matrixCO2['SOLAR_ST'],err_matrixCO2['T_ST'],err_matrixCO2['SPEC_ST'],err_matrixCO2['NOISE_ST'],err_matrixCO2['BASE_SY'],err_matrixCO2['ILS_SY'],err_matrixCO2['LOS_SY'],err_matrixCO2['SOLAR_SY'],err_matrixCO2['T_SY'],err_matrixCO2['SPEC_SY'],err_matrixCO2['NOISE_SY'])
			print err_matrixCO2['BASE_ST'],err_matrixCO2['ILS_ST'],err_matrixCO2['LOS_ST'],err_matrixCO2['SOLAR_ST'],err_matrixCO2['T_ST'],err_matrixCO2['SPEC_ST'],err_matrixCO2['NOISE_ST'],err_matrixCO2['BASE_SY'],err_matrixCO2['ILS_SY'],err_matrixCO2['LOS_SY'],err_matrixCO2['SOLAR_SY'],err_matrixCO2['T_SY'],err_matrixCO2['SPEC_SY'],err_matrixCO2['NOISE_SY']
			#print 'CO2 TOTAL ST = %.06f, SY = %.06f' % (np.sqrt(err_matrixCO2['BASE_ST']**2 + err_matrixCO2['ILS_ST']**2 + err_matrixCO2['LOS_ST']**2 + err_matrixCO2['SOLAR_ST']**2 + err_matrixCO2['T_ST']**2 + err_matrixCO2['SPEC_ST']**2 + err_matrixCO2['NOISE_ST']**2), np.sqrt(err_matrixCO2['BASE_SY']**2 + err_matrixCO2['ILS_SY']**2 + err_matrixCO2['LOS_SY']**2 + err_matrixCO2['SOLAR_SY']**2 + err_matrixCO2['T_SY']**2 + err_matrixCO2['SPEC_SY']**2 + err_matrixCO2['NOISE_SY']**2))
			print np.sqrt(err_matrixCO2['BASE_ST']**2 + err_matrixCO2['ILS_ST']**2 + err_matrixCO2['LOS_ST']**2 + err_matrixCO2['SOLAR_ST']**2 + err_matrixCO2['T_ST']**2 + err_matrixCO2['SPEC_ST']**2 + err_matrixCO2['NOISE_ST']**2), np.sqrt(err_matrixCO2['BASE_SY']**2 + err_matrixCO2['ILS_SY']**2 + err_matrixCO2['LOS_SY']**2 + err_matrixCO2['SOLAR_SY']**2 + err_matrixCO2['T_SY']**2 + err_matrixCO2['SPEC_SY']**2 + err_matrixCO2['NOISE_SY']**2)
			#print 'O2: BASE ST = %.06f, ILS ST = %.06e, LOS ST = %.06f, SOLAR ST = %.06f, T ST = %.06f, SPEC ST = %.06f, NOISE ST = %.06f, BASE SY = %.06f, ILS SY = %.06e, LOS SY = %.06f, SOLAR SY = %.06f, T SY = %.06f, SPEC SY = %.06f, NOISE SY = %.06f' % (err_matrixO2['BASE_ST'],err_matrixO2['ILS_ST'],err_matrixO2['LOS_ST'],err_matrixO2['SOLAR_ST'],err_matrixO2['T_ST'],err_matrixO2['SPEC_ST'],err_matrixO2['NOISE_ST'],err_matrixO2['BASE_SY'],err_matrixO2['ILS_SY'],err_matrixO2['LOS_SY'],err_matrixO2['SOLAR_SY'],err_matrixO2['T_SY'],err_matrixO2['SPEC_SY'],err_matrixO2['NOISE_SY'])
			print err_matrixO2['BASE_ST'],err_matrixO2['ILS_ST'],err_matrixO2['LOS_ST'],err_matrixO2['SOLAR_ST'],err_matrixO2['T_ST'],err_matrixO2['SPEC_ST'],err_matrixO2['NOISE_ST'],err_matrixO2['BASE_SY'],err_matrixO2['ILS_SY'],err_matrixO2['LOS_SY'],err_matrixO2['SOLAR_SY'],err_matrixO2['T_SY'],err_matrixO2['SPEC_SY'],err_matrixO2['NOISE_SY']
			#print 'O2 TOTAL ST = %.06f, SY = %.06f' % (np.sqrt(err_matrixO2['BASE_ST']**2 + err_matrixO2['ILS_ST']**2 + err_matrixO2['LOS_ST']**2 + err_matrixO2['SOLAR_ST']**2 + err_matrixO2['T_ST']**2 + err_matrixO2['SPEC_ST']**2 + err_matrixO2['NOISE_ST']**2), np.sqrt(err_matrixO2['BASE_SY']**2 + err_matrixO2['ILS_SY']**2 + err_matrixO2['LOS_SY']**2 + err_matrixO2['SOLAR_SY']**2 + err_matrixO2['T_SY']**2 + err_matrixO2['SPEC_SY']**2 + err_matrixO2['NOISE_SY']**2))
			print np.sqrt(err_matrixO2['BASE_ST']**2 + err_matrixO2['ILS_ST']**2 + err_matrixO2['LOS_ST']**2 + err_matrixO2['SOLAR_ST']**2 + err_matrixO2['T_ST']**2 + err_matrixO2['SPEC_ST']**2 + err_matrixO2['NOISE_ST']**2), np.sqrt(err_matrixO2['BASE_SY']**2 + err_matrixO2['ILS_SY']**2 + err_matrixO2['LOS_SY']**2 + err_matrixO2['SOLAR_SY']**2 + err_matrixO2['T_SY']**2 + err_matrixO2['SPEC_SY']**2 + err_matrixO2['NOISE_SY']**2)
			#print 'POINTING ERROR '+bms[ii]+' (SZA = %.06f) = %.06f' % (err_matrixCO2['SZA'],np.sqrt((err_matrixCO2['LOS_ST']+err_matrixCO2['LOS_SY'])**2 + (err_matrixO2['LOS_ST']+err_matrixO2['LOS_SY'])**2))
			#print 'XCO2 NOISE '+bms[ii]+' (SZA = %.06f) = %.06f' % (err_matrixCO2['SZA'],np.sqrt(err_matrixCO2['NOISE_ST']**2 + err_matrixO2['NOISE_ST']**2))

			'''
			plt.figure(1)
			plt.suptitle('CO2 STATISTICAL ERRORS')
			plt.yscale('log')
			plt.plot(tm_CO2,err_matrixCO2['BASE_ST'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(tm_CO2,err_matrixCO2['ILS_ST'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['LOS_ST'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['SOLAR_ST'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['T_ST'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			#plt.plot(tm_CO2,err_matrixCO2['SPEC_ST'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['NOISE_ST'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')		
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('CO2_ST.eps')
		
			plt.figure(2)
			plt.suptitle('CO2 SYSTEMATIC ERRORS')
			plt.yscale('log')
			plt.plot(tm_CO2,err_matrixCO2['BASE_SY'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(tm_CO2,err_matrixCO2['ILS_SY'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['LOS_SY'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['SOLAR_SY'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['T_SY'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			plt.plot(tm_CO2,err_matrixCO2['SPEC_SY'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			#plt.plot(tm_CO2,err_matrixCO2['NOISE_SY'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('CO2_SY.eps')

			plt.figure(3)
			plt.suptitle('O2 STATISTICAL ERRORS')
			plt.yscale('log')
			plt.plot(tm_O2,err_matrixO2['BASE_ST'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(tm_O2,err_matrixO2['ILS_ST'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['LOS_ST'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['SOLAR_ST'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['T_ST'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			#plt.plot(tm_O2,err_matrixO2['SPEC_ST'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['NOISE_ST'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('O2_ST.eps')

			plt.figure(4)
			plt.suptitle('O2 SYSTEMATIC ERRORS')
			plt.yscale('log')
			plt.plot(tm_O2,err_matrixO2['BASE_SY'],label='BASE '+bms[ii],linestyle=lsty[ii],color='b',linewidth=lw_n)
			#plt.plot(tm_O2,err_matrixO2['ILS_SY'],label='ILS '+bms[ii],linestyle=lsty[ii],color='c',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['LOS_SY'],label='LOS '+bms[ii],linestyle=lsty[ii],color='g',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['SOLAR_SY'],label='SOLAR '+bms[ii],linestyle=lsty[ii],color='y',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['T_SY'],label='T '+bms[ii],linestyle=lsty[ii],color='k',linewidth=lw_n)
			plt.plot(tm_O2,err_matrixO2['SPEC_SY'],label='SPEC '+bms[ii],linestyle=lsty[ii],color='m',linewidth=lw_n)
			#plt.plot(tm_O2,err_matrixO2['NOISE_SY'],label='NOISE '+bms[ii],linestyle=lsty[ii],color='r',linewidth=lw_n)
			plt.xlabel('Local Time')
			plt.ylabel('%')
			#plt.xticks(rotation=45)
			plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
			plt.legend(ncol=2,loc=8,prop={'size':8},handlelength=3)
			#plt.savefig('O2_SY.eps')
			'''

	else:
		print 'PATHS:',testCO2,testO2
		err_matrixCO2 = lee(testCO2,fec,'CO2')	### LLAMA A LEE()
		err_matrixO2 = lee(testO2,fec,'O2')	### LLAMA A LEE()

		print err_matrixCO2[0]
		np.savetxt(fec+'_CO2_err.txt',err_matrixCO2,fmt="%s")#"%s %i %.4f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e")
		np.savetxt(fec+'_O2_err.txt',err_matrixO2,fmt="%s")#"%s %i %.4f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e")


plt.show()
