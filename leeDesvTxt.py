import numpy as np
import sys
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates
import os
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'


# LEE ARCHIVOS DESV_****.TXT Y GRAFICA VALORES

##############################

def julian_day(fecha):
    return 367*fecha.year-(7*(fecha.year+((fecha.month+9)/12))/4)+(275*fecha.month/9)+fecha.day+1721013.5-0.5*np.sign((100*fecha.year)+fecha.month-190002.5)+0.5+fecha.hour/24.0+fecha.minute/(60.0*24.0)+fecha.second/(3600.0*24.0)

##############################

def sun_angle(jt_UT):
	import ephem
	lat = 19.1187
	lon = -98.6552
	alt = 3985.0
	o = ephem.Observer()
	o.lat = lat*ephem.degree
	o.long =lon*ephem.degree
	o.elev=alt
	o.date = ephem.Date(jt_UT-2415020)
	sun = ephem.Sun(o)
       
	elev=float(sun.alt)
	azimuth=float(sun.az)
        astro_SZA = np.abs(90.0-(elev*180/np.pi))
        astro_SAA = azimuth*180/np.pi
        #dn=0.000292*np.exp(-alt/6000.0) # pressure actualisation
        #SZA_refacted = 180.0/np.pi*np.arcsin(np.sin(np.pi*astro_SZA/180.0)/(1.0+dn))
	return astro_SZA

##############################

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

##############################

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

##############################

def DESV(especie):

	f_name = 'DESV_'+especie+'.txt'
	f_txt = open(f_name,'r')
	f_lines = f_txt.readlines()
	f_txt.close()

	f_fechasCaF2 = np.genfromtxt('fechasCaF2.txt',dtype=str)

	mat_desv = []
	for ii,line in enumerate(f_lines):
		if line[0:5] == 'FECHA':
			line_cachos = line.strip().split()
			numDeMed = int(line_cachos[5])
			fecha = line_cachos[1]
			resto = numDeMed
			jj = 1
			while resto != 0:
				firstlinesplit = f_lines[ii+jj].strip().split()
				secondlinesplit = f_lines[ii+jj+1].strip().split()
				medCont = int(firstlinesplit[2].split(')')[0][1:])			### NUMERO DE MEDICIONES CONTIGUAS
				resto = resto - medCont
				primerahora = dt.datetime.strptime(fecha+firstlinesplit[3],'%y%m%d%H:%M:%S')	### TIEMPO DE LA PRIMERA MEDICION USADA
				ultimahora = dt.datetime.strptime(fecha+firstlinesplit[-1],'%y%m%d%H:%M:%S')	### TIEMPO DE LA ULTIMA MEDICION USADA
				deltahora = (ultimahora-primerahora).total_seconds()			### DURACION DE MEDICIONES EN SEGUNDOS
				meanhora = primerahora + dt.timedelta(seconds = deltahora/2,hours = +6)	### HORA PARA CALCULAR EL SZA DE REFERENCIA EN UT
				meanjd = julian_day(meanhora)						### JULIAN DAY
				meansza = sun_angle(meanjd)						### ASTRONOMICAL SOLAR ZENITH ANGLE
				cont_mean = float(secondlinesplit[1])						### PROMEDIO DE MEDICIONES
				cont_std = float(secondlinesplit[3])						### DESVIACION DE MEDICIONES
				cont_percent = float(secondlinesplit[5])					### DESVIACION/PROMEDIO
				bms = 'kbr'
				fcaf2_con = (fecha == f_fechasCaF2[0])
				for fec in f_fechasCaF2[1:]:
					fcaf2_con = fcaf2_con | (fecha == fec)
				if (fcaf2_con == True):
					bms = 'caf2'
				mat_desv.append([fecha,numDeMed,medCont,primerahora,ultimahora,deltahora,meansza,cont_mean,cont_std,cont_percent,bms])
				jj = jj + 2
			jj = 1

	mat_desv = np.array(mat_desv)
	dt_desv = mat_desv[:,3]
	my_desv = np.array(['%04i-%02i' % (ele.year,ele.month) for ele in dt_desv])
	kbr_con = ((mat_desv[:,8] != 0) & (mat_desv[:,10] == 'kbr'))
	caf2_con = ((mat_desv[:,8] != 0) & (mat_desv[:,10] == 'caf2'))
	kbr02_con = ((mat_desv[:,8] != 0) & (mat_desv[:,2] <= 10) & (mat_desv[:,10] == 'kbr') & (my_desv == '2014-02'))
	caf202_con = ((mat_desv[:,8] != 0) & (mat_desv[:,2] <= 10) & (mat_desv[:,10] == 'caf2') & (my_desv == '2014-02'))

	print especie
	print 'Total:',np.mean(mat_desv[:,9]),len(mat_desv[:,9]),'( STD == 0:',len(mat_desv[mat_desv[:,8] == 0][:,9]),')'
	print 'KBr:',np.mean(mat_desv[kbr_con][:,9]),len(mat_desv[kbr_con][:,9]),'MEAS:',np.sum(mat_desv[kbr_con][:,2])
	print 'CaF2:',np.mean(mat_desv[caf2_con][:,9]),len(mat_desv[caf2_con][:,9]),'MEAS:',np.sum(mat_desv[caf2_con][:,2])
	print '14-02:',np.mean(mat_desv[my_desv == '2014-02'][:,9]),len(mat_desv[my_desv == '2014-02'][:,9]),'( # > 10 | STD == 0:',len(mat_desv[((my_desv == '2014-02') & (mat_desv[:,2] > 10)) | ((my_desv == '2014-02') & (mat_desv[:,8] == 0))] [:,9]),')'
	print 'KBr (14-02):',np.mean(mat_desv[kbr02_con][:,9]),len(mat_desv[kbr02_con][:,9])
	print 'CaF2 (14-02):',np.mean(mat_desv[caf202_con][:,9]),len(mat_desv[caf202_con][:,9])

	mc_hist, mc_be = np.histogram(mat_desv[:,2],bins=50)					### NUMERO DE MEDICIONES CONTIGUAS
	mc_bp = np.array([(ele + mc_be[ii+1])/2 for ii,ele in enumerate(mc_be[:-1])])

	dhr_hist,dhr_be = np.histogram(mat_desv[:,5],bins=50)					### DELTA HORA
	dhr_bp = np.array([(ele + dhr_be[ii+1])/2 for ii,ele in enumerate(dhr_be[:-1])])

	cs_hist, cs_be = np.histogram(mat_desv[mat_desv[:,8] != 0][:,8],bins=50)		### DESVIACION ESTANDAR
	cs_bp = np.array([(ele + cs_be[ii+1])/2 for ii,ele in enumerate(cs_be[:-1])])

	if especie == 'XCO2':
		cp_hist, cp_be = np.histogram(mat_desv[mat_desv[:,8] != 0][:,9],bins=25)		### PERCENT (SIN STD = 0)
	else:
		cp_hist, cp_be = np.histogram(mat_desv[mat_desv[:,8] != 0][:,9],bins=100)		### PERCENT (SIN STD = 0)
	cp_bp = np.array([(ele + cp_be[ii+1])/2 for ii,ele in enumerate(cp_be[:-1])])
	cp_hist = np.array([float(ele)/float(np.sum(cp_hist)) for ele in cp_hist])

	cp_kbr_hist, cp_kbr_be = np.histogram(mat_desv[kbr_con][:,9],bins=cp_be)		### PERCENT KBr (SIN STD = 0)
	#cp_kbr_bp = np.array([(ele + cp_kbr_be[ii+1])/2 for ii,ele in enumerate(cp_kbr_be[:-1])])
	cp_kbr_hist = np.array([float(ele)/float(np.sum(cp_kbr_hist)) for ele in cp_kbr_hist])

	cp_caf2_hist, cp_caf2_be = np.histogram(mat_desv[caf2_con][:,9],bins=cp_be)	### PERCENT CaF2 (SIN STD = 0)
	#cp_caf2_bp = np.array([(ele + cp_caf2_be[ii+1])/2 for ii,ele in enumerate(cp_caf2_be[:-1])])
	cp_caf2_hist = np.array([float(ele)/float(np.sum(cp_caf2_hist)) for ele in cp_caf2_hist])

	if especie == 'XCO2':
		cp02_hist, cp02_be = np.histogram(mat_desv[(mat_desv[:,8] != 0) & (my_desv == '2014-02')][:,9],bins=10)	### PERCENT 14-02 (SIN STD = 0)
	else:
		cp02_hist, cp02_be = np.histogram(mat_desv[(mat_desv[:,8] != 0) & (my_desv == '2014-02')][:,9],bins=25)	### PERCENT 14-02 (SIN STD = 0)
	cp02_bp = np.array([(ele + cp02_be[ii+1])/2 for ii,ele in enumerate(cp02_be[:-1])])
	cp02_hist = np.array([float(ele)/float(np.sum(cp02_hist)) for ele in cp02_hist])

	cp_kbr02_hist, cp_kbr02_be = np.histogram(mat_desv[kbr02_con][:,9],bins=cp02_be)		### PERCENT KBr 14-02 (SIN STD = 0)
	#cp_kbr02_bp = np.array([(ele + cp_kbr02_be[ii+1])/2 for ii,ele in enumerate(cp_kbr02_be[:-1])])
	cp_kbr02_hist = np.array([float(ele)/float(np.sum(cp_kbr02_hist)) for ele in cp_kbr02_hist])

	cp_caf202_hist, cp_caf202_be = np.histogram(mat_desv[caf202_con][:,9],bins=cp02_be)	### PERCENT CaF2 14-02 (SIN STD = 0)
	#cp_caf202_bp = np.array([(ele + cp_caf202_be[ii+1])/2 for ii,ele in enumerate(cp_caf202_be[:-1])])
	cp_caf202_hist = np.array([float(ele)/float(np.sum(cp_caf202_hist)) for ele in cp_caf202_hist])

	plt.figure()	
	plt.suptitle(especie+' TODOS')
	#plt.subplot(311)
	plt.plot(cp_bp,cp_hist,'-k',label=None)
	#plt.xlim(0,1)
	#plt.subplot(312)
	plt.plot(cp_bp,cp_kbr_hist,'-r',label='KBr')
	#plt.xlim(0,1)
	#plt.subplot(313)
	plt.plot(cp_bp,cp_caf2_hist,'-b',label='CaF2')
	plt.xlabel('%')
	if especie != 'XCO2':
		plt.xlim(0,0.5)
	else:
		plt.xlim(0,0.3)
	plt.grid()
	plt.legend()

	plt.figure()	
	plt.suptitle(especie+' 14-02')
	#plt.subplot(311)
	plt.plot(cp02_bp,cp02_hist,'-k',label=None)
	#plt.xlim(0,0.4)
	#plt.subplot(312)
	plt.plot(cp02_bp,cp_kbr02_hist,'-r',label='KBr')
	#plt.xlim(0,0.4)
	#plt.subplot(313)
	plt.plot(cp02_bp,cp_caf202_hist,'-b',label='CaF2')
	plt.xlabel('%')
	if especie != 'XCO2':
		plt.xlim(0,0.5)
	else:
		plt.xlim(0,0.3)
	plt.grid()
	plt.legend()

	plt.figure()
	plt.suptitle(especie+' TODOS')
	plt.plot(mat_desv[kbr_con][:,6],mat_desv[kbr_con][:,9],'+r',label='KBr')
	plt.plot(mat_desv[caf2_con][:,6],mat_desv[caf2_con][:,9],'xb',label='CaF2')
	plt.xlabel('SZA')
	plt.ylabel('%')
	plt.legend()

	plt.figure()
	plt.suptitle(especie+' 14-02')
	plt.plot(mat_desv[kbr02_con][:,6],mat_desv[kbr02_con][:,9],'+r',label='KBr')
	plt.plot(mat_desv[caf202_con][:,6],mat_desv[caf202_con][:,9],'xb',label='CaF2')
	plt.xlabel('SZA')
	plt.ylabel('%')
	plt.legend()

##############################

def COLUMNS(archivo_txt):

	def HISTOS_COL(archivo_txt,f_col,KBr_idx,CaF2_idx,numbins=25):

		hist_CO2_KBr,be_CO2_KBr = np.histogram(f_col['totcolCO2'][KBr_idx],numbins)
		hist_CO2_norm_KBr = np.array([float(ele)/float(np.sum(hist_CO2_KBr)) for ele in hist_CO2_KBr])
		bp_CO2_KBr = np.array([(ele + be_CO2_KBr[ii+1])/2 for ii,ele in enumerate(be_CO2_KBr[:-1])])
		hist_CO2_CaF2,be_CO2_CaF2 = np.histogram(f_col['totcolCO2'][CaF2_idx],be_CO2_KBr)
		hist_CO2_norm_CaF2 = np.array([float(ele)/float(np.sum(hist_CO2_CaF2)) for ele in hist_CO2_CaF2])
		bp_CO2_CaF2 = np.array([(ele + be_CO2_CaF2[ii+1])/2 for ii,ele in enumerate(be_CO2_CaF2[:-1])])

		dtr_hist_CO2_KBr,dtr_be_CO2_KBr = np.histogram(dtr_CO2[KBr_idx],numbins)
		dtr_hist_CO2_norm_KBr = np.array([float(ele)/float(np.sum(dtr_hist_CO2_KBr)) for ele in dtr_hist_CO2_KBr])
		dtr_bp_CO2_KBr = np.array([(ele + dtr_be_CO2_KBr[ii+1])/2 for ii,ele in enumerate(dtr_be_CO2_KBr[:-1])])
		dtr_hist_CO2_CaF2,dtr_be_CO2_CaF2 = np.histogram(dtr_CO2[CaF2_idx],dtr_be_CO2_KBr)
		dtr_hist_CO2_norm_CaF2 = np.array([float(ele)/float(np.sum(dtr_hist_CO2_CaF2)) for ele in dtr_hist_CO2_CaF2])
		dtr_bp_CO2_CaF2 = np.array([(ele + dtr_be_CO2_CaF2[ii+1])/2 for ii,ele in enumerate(dtr_be_CO2_CaF2[:-1])])

		hist_O2_KBr,be_O2_KBr = np.histogram(f_col['totcolO2'][KBr_idx],numbins)
		hist_O2_norm_KBr = np.array([float(ele)/float(np.sum(hist_O2_KBr)) for ele in hist_O2_KBr])
		bp_O2_KBr = np.array([(ele + be_O2_KBr[ii+1])/2 for ii,ele in enumerate(be_O2_KBr[:-1])])
		hist_O2_CaF2,be_O2_CaF2 = np.histogram(f_col['totcolO2'][CaF2_idx],be_O2_KBr)
		hist_O2_norm_CaF2 = np.array([float(ele)/float(np.sum(hist_O2_CaF2)) for ele in hist_O2_CaF2])
		bp_O2_CaF2 = np.array([(ele + be_O2_CaF2[ii+1])/2 for ii,ele in enumerate(be_O2_CaF2[:-1])])

		dtr_hist_O2_KBr,dtr_be_O2_KBr = np.histogram(dtr_O2[KBr_idx],numbins)
		dtr_hist_O2_norm_KBr = np.array([float(ele)/float(np.sum(dtr_hist_O2_KBr)) for ele in dtr_hist_O2_KBr])
		dtr_bp_O2_KBr = np.array([(ele + dtr_be_O2_KBr[ii+1])/2 for ii,ele in enumerate(dtr_be_O2_KBr[:-1])])
		dtr_hist_O2_CaF2,dtr_be_O2_CaF2 = np.histogram(dtr_O2[CaF2_idx],dtr_be_O2_KBr)
		dtr_hist_O2_norm_CaF2 = np.array([float(ele)/float(np.sum(dtr_hist_O2_CaF2)) for ele in dtr_hist_O2_CaF2])
		dtr_bp_O2_CaF2 = np.array([(ele + dtr_be_O2_CaF2[ii+1])/2 for ii,ele in enumerate(dtr_be_O2_CaF2[:-1])])

		hist_XCO2_KBr,be_XCO2_KBr = np.histogram(f_col['XCO2'][KBr_idx],numbins)
		hist_XCO2_norm_KBr = np.array([float(ele)/float(np.sum(hist_XCO2_KBr)) for ele in hist_XCO2_KBr])
		bp_XCO2_KBr = np.array([(ele + be_XCO2_KBr[ii+1])/2 for ii,ele in enumerate(be_XCO2_KBr[:-1])])
		hist_XCO2_CaF2,be_XCO2_CaF2 = np.histogram(f_col['XCO2'][CaF2_idx],be_XCO2_KBr)
		hist_XCO2_norm_CaF2 = np.array([float(ele)/float(np.sum(hist_XCO2_CaF2)) for ele in hist_XCO2_CaF2])
		bp_XCO2_CaF2 = np.array([(ele + be_XCO2_CaF2[ii+1])/2 for ii,ele in enumerate(be_XCO2_CaF2[:-1])])

		dtr_hist_XCO2_KBr,dtr_be_XCO2_KBr = np.histogram(dtr_XCO2[KBr_idx],numbins)
		dtr_hist_XCO2_norm_KBr = np.array([float(ele)/float(np.sum(dtr_hist_XCO2_KBr)) for ele in dtr_hist_XCO2_KBr])
		dtr_bp_XCO2_KBr = np.array([(ele + dtr_be_XCO2_KBr[ii+1])/2 for ii,ele in enumerate(dtr_be_XCO2_KBr[:-1])])
		dtr_hist_XCO2_CaF2,dtr_be_XCO2_CaF2 = np.histogram(dtr_XCO2[CaF2_idx],dtr_be_XCO2_KBr)
		dtr_hist_XCO2_norm_CaF2 = np.array([float(ele)/float(np.sum(dtr_hist_XCO2_CaF2)) for ele in dtr_hist_XCO2_CaF2])
		dtr_bp_XCO2_CaF2 = np.array([(ele + dtr_be_XCO2_CaF2[ii+1])/2 for ii,ele in enumerate(dtr_be_XCO2_CaF2[:-1])])

		plt.figure()
		plt.suptitle('CO2 ('+archivo_txt+')')
		plt.subplot(211)
		plt.plot(bp_CO2_CaF2,hist_CO2_norm_CaF2,'-b',label='CaF2')
		plt.plot(bp_CO2_KBr,hist_CO2_norm_KBr,'-r',label='KBr')
		plt.xlabel('CO2 (molec/cm2)')
		#plt.ylim(0.0,0.35)
		plt.legend()
		plt.subplot(212)
		plt.plot(dtr_bp_CO2_CaF2,dtr_hist_CO2_norm_CaF2,'-b',label='CaF2')
		plt.plot(dtr_bp_CO2_KBr,dtr_hist_CO2_norm_KBr,'-r',label='KBr')
		plt.xlabel('DETRENDED CO2 (molec/cm2)')
		#plt.ylim(0.0,0.35)
		plt.legend()

		plt.figure()
		plt.suptitle('O2 ('+archivo_txt+')')
		plt.subplot(211)
		plt.plot(bp_O2_CaF2,hist_O2_norm_CaF2,'-b',label='CaF2')
		plt.plot(bp_O2_KBr,hist_O2_norm_KBr,'-r',label='KBr')
		plt.xlabel('O2 (molec/cm2)')
		#plt.ylim(0.0,0.35)
		plt.legend()
		plt.subplot(212)
		plt.plot(dtr_bp_O2_CaF2,dtr_hist_O2_norm_CaF2,'-b',label='CaF2')
		plt.plot(dtr_bp_O2_KBr,dtr_hist_O2_norm_KBr,'-r',label='KBr')
		plt.xlabel('DETRENDED O2 (molec/cm2)')
		#plt.ylim(0.0,0.35)
		plt.legend()

		plt.figure()
		plt.suptitle('XCO2 ('+archivo_txt+')')
		plt.subplot(211)
		plt.plot(bp_XCO2_CaF2,hist_XCO2_norm_CaF2,'-b',label='CaF2')
		plt.plot(bp_XCO2_KBr,hist_XCO2_norm_KBr,'-r',label='KBr')
		plt.xlabel('XCO2 (ppm)')
		#plt.ylim(0.0,0.35)
		plt.legend()
		plt.subplot(212)
		plt.plot(dtr_bp_XCO2_CaF2,dtr_hist_XCO2_norm_CaF2,'-b',label='CaF2')
		plt.plot(dtr_bp_XCO2_KBr,dtr_hist_XCO2_norm_KBr,'-r',label='KBr')
		plt.xlabel('DETRENDED XCO2 (ppm)')
		#plt.ylim(0.0,0.35)
		plt.legend()

	##############################

	def ENSAMBLE(archivo_txt,f_col,KBr_idx,CaF2_idx,flag_todo=0,modo=0):	### MODO: 1: FECHAS SIMILARES EN DIFERENTES ANOS, 0: FECHAS MAS RECIENTES
		### flag_todo: 0 PARA USAR UN ENSAMBLE IGUAL AL DE CaF2, 1 PARA USAR TODOS LOS DATOS DE KBr
	
		dt_txt = np.dtype([('fname',np.str_,30),('epochTime',int),('SZA',float), \
		('XCO2',float),('totcolCO2',float),('totcolO2',float),('LOSerrCO2',float),('NOISEerrCO2',float), \
		('LOSerrO2',float),('NOISEerrO2',float)])
		
		KBr_en = f_col[KBr_idx]
		KBr_dt = np.array([dt.datetime.fromtimestamp(ele) for ele in KBr_en['epochTime']])

		CaF2_en = f_col[CaF2_idx]
		CaF2_yr = np.array([dt.datetime.fromtimestamp(ele).year for ele in CaF2_en['epochTime']])
		### DIFERENCIA DE SEGUNDOS ENTRE FECHAS DE KBr Y FECHAS DE REFERENCIA DE 2014 Y 2015 DE CaF2
		KBr_ts_2014 = np.array([abs((ele.replace(year=2001) - dt.datetime(2001,02,18)).total_seconds()) for ele in KBr_dt])	### SE USO 18-02-2014 DE REFERENCIA
		KBr_ts_2015 = np.array([abs((ele.replace(year=2001) - dt.datetime(2001,05,18)).total_seconds()) for ele in KBr_dt])	### SE USO 18-05-2015 DE REFERENCIA
		hist_SZA_KBr, be_SZA_KBr = np.histogram(KBr_en['SZA'],160)
		hist_SZA_norm_KBr = np.array([float(ele)/float(np.sum(hist_SZA_KBr)) for ele in hist_SZA_KBr])
		bp_SZA_KBr = np.array([(ele + be_SZA_KBr[ii+1])/2 for ii,ele in enumerate(be_SZA_KBr[:-1])])

		hist_SZA_CaF2, be_SZA_CaF2 = np.histogram(CaF2_en['SZA'],be_SZA_KBr)
		hist_SZA_norm_CaF2 = np.array([float(ele)/float(np.sum(hist_SZA_CaF2)) for ele in hist_SZA_CaF2])
		bp_SZA_CaF2 = np.array([(ele + be_SZA_CaF2[ii+1])/2 for ii,ele in enumerate(be_SZA_CaF2[:-1])])

		### DOS HISTOGRAMAS, UNO PARA 2014 Y OTRO PARA 2015
		bool_exp_2014 = (CaF2_yr == 2014)
		bool_exp_2015 = (CaF2_yr == 2015)

		print 'CaF2 2014: %i, 2015: %i' % (len(CaF2_en[bool_exp_2014]),len(CaF2_en[bool_exp_2015]))

		hist_SZA_CaF2_2014, be_SZA_CaF2_2014 = np.histogram(CaF2_en[bool_exp_2014]['SZA'],be_SZA_KBr)
		hist_SZA_norm_CaF2_2014 = np.array([float(ele)/float(np.sum(hist_SZA_CaF2_2014)) for ele in hist_SZA_CaF2_2014])
		bp_SZA_CaF2_2014 = np.array([(ele + be_SZA_CaF2_2014[ii+1])/2 for ii,ele in enumerate(be_SZA_CaF2_2014[:-1])])
		
		hist_SZA_CaF2_2015, be_SZA_CaF2_2015 = np.histogram(CaF2_en[bool_exp_2015]['SZA'],be_SZA_KBr)
		hist_SZA_norm_CaF2_2015 = np.array([float(ele)/float(np.sum(hist_SZA_CaF2_2015)) for ele in hist_SZA_CaF2_2015])
		bp_SZA_CaF2_2015 = np.array([(ele + be_SZA_CaF2_2015[ii+1])/2 for ii,ele in enumerate(be_SZA_CaF2_2015[:-1])])

		print 'LEN ENSAMBLE KBr',len(KBr_en),np.sum(hist_SZA_KBr)
		print 'LEN ENSAMBLE CaF2',len(CaF2_en),np.sum(hist_SZA_CaF2_2014)+np.sum(hist_SZA_CaF2_2015)

		if ((os.path.exists('CaF2ensamble.txt')) & (os.path.exists('KBrensamble%i.txt' % (flag_todo)))):

			print 'ARCHIVOS CaF2ensamble.txt y KBrensamble%i.txt ENCONTRADOS :)' % (flag_todo)
			CaF2_ensamble = np.genfromtxt('CaF2ensamble.txt',dtype=dt_txt,delimiter = " ")
			KBr_ensamble = np.genfromtxt('KBrensamble%i.txt' % flag_todo,dtype=dt_txt,delimiter = " ")

		else:
			
			print 'ARCHIVOS txt NO ENCONTRADOS :('

			CaF2_ensamble = CaF2_en		### ENSAMBLE DE CaF2 CON VALORES DE SZA SIMILARES A CAF2
			if flag_todo == 0:		
				KBr_ensamble = np.empty([np.sum(hist_SZA_CaF2_2014)+np.sum(hist_SZA_CaF2_2015)],dtype=dt_txt)		### ENSAMBLE DE KBr CON VALORES DE SZA SIMILARES A CAF2

				if modo == 1:
					dum_CaF2 = 0
					for idx,bedges in enumerate(be_SZA_CaF2_2014[:-1]):
						exp_bool_SZA_CaF2 = ((CaF2_ensamble['SZA'] >= bedges) & (CaF2_ensamble['SZA'] <= be_SZA_CaF2_2014[idx+1]))
						exp_bool_SZA_KBr = ((KBr_en['SZA'] >= bedges) & (KBr_en['SZA'] <= be_SZA_CaF2_2014[idx+1]))
						ln_CaF2 = hist_SZA_CaF2_2014[idx]
						sm_KBr_ts_2014 = KBr_ts_2014[exp_bool_SZA_KBr]
						KBr_ts_idx_2014 = np.argsort(sm_KBr_ts_2014)
						sm_KBr_ensamble = KBr_en[exp_bool_SZA_KBr]
						sm_KBr_ensamble = sm_KBr_ensamble[KBr_ts_idx_2014]

						print 'IN 2014 FOR %.4f <= SZA <= %.4f, KBr SIZE: %i, CaF2 SIZE: %i' % (bedges,be_SZA_CaF2_2014[idx+1],len(sm_KBr_ensamble),ln_CaF2)
						#print >> open('KBr_ALL_2014.txt','a'),'IN 2014 FOR %.4f <= SZA <= %.4f, KBr SIZE: %i, CaF2 SIZE: %i' % (bedges,be_SZA_CaF2_2014[idx+1],len(sm_KBr_ensamble),ln_CaF2),sm_KBr_ensamble['fname']
					
						for rr in range(0,ln_CaF2):
							KBr_ensamble[dum_CaF2+rr-1] = sm_KBr_ensamble[rr]
						dum_CaF2 = dum_CaF2 + ln_CaF2

					for idx,bedges in enumerate(be_SZA_CaF2_2015[:-1]):
						exp_bool_SZA_CaF2 = ((CaF2_ensamble['SZA'] >= bedges) & (CaF2_ensamble['SZA'] <= be_SZA_CaF2_2015[idx+1]))
						exp_bool_SZA_KBr = ((KBr_en['SZA'] >= bedges) & (KBr_en['SZA'] <= be_SZA_CaF2_2015[idx+1]))
						ln_CaF2 = hist_SZA_CaF2_2015[idx]
						sm_KBr_ts_2015 = KBr_ts_2015[exp_bool_SZA_KBr]
						KBr_ts_idx_2015 = np.argsort(sm_KBr_ts_2015)
						sm_KBr_ensamble = KBr_en[exp_bool_SZA_KBr]
						sm_KBr_ensamble = sm_KBr_ensamble[KBr_ts_idx_2015]

						print 'IN 2015 FOR %.4f <= SZA <= %.4f, KBr SIZE: %i, CaF2 SIZE: %i' % (bedges,be_SZA_CaF2_2015[idx+1],len(sm_KBr_ensamble),ln_CaF2)
						print >> open('KBr_ALL_2015.txt','a'),'IN 2015 FOR %.4f <= SZA <= %.4f, KBr SIZE: %i, CaF2 SIZE: %i' % (bedges,be_SZA_CaF2_2015[idx+1],len(sm_KBr_ensamble),ln_CaF2),sm_KBr_ts_2015#sm_KBr_ensamble['fname']
					
						for rr in range(0,ln_CaF2):
							KBr_ensamble[dum_CaF2+rr-1] = sm_KBr_ensamble[rr]
						dum_CaF2 = dum_CaF2 + ln_CaF2
				elif modo == 0:
					dum_CaF2 = 0
					for idx,bedges in enumerate(be_SZA_CaF2[:-1]):
						exp_bool_SZA_CaF2 = ((CaF2_ensamble['SZA'] >= bedges) & (CaF2_ensamble['SZA'] <= be_SZA_CaF2[idx+1]))
						exp_bool_SZA_KBr = ((KBr_en['SZA'] >= bedges) & (KBr_en['SZA'] <= be_SZA_CaF2[idx+1]))
						ln_CaF2 = hist_SZA_CaF2[idx]
						sm_KBr_ensamble = KBr_en[exp_bool_SZA_KBr]
						sm_KBr_dt = KBr_dt[exp_bool_SZA_KBr]
						idx_KBr_dt = np.argsort(sm_KBr_dt)
						sm_KBr_ensamble = np.flipud(sm_KBr_ensamble[idx_KBr_dt])

						print 'FOR %.4f <= SZA <= %.4f, KBr SIZE: %i, CaF2 SIZE: %i' % (bedges,be_SZA_CaF2[idx+1],len(sm_KBr_ensamble),ln_CaF2)		
						for rr in range(0,ln_CaF2):
							KBr_ensamble[dum_CaF2+rr-1] = sm_KBr_ensamble[rr]
						dum_CaF2 = dum_CaF2 + ln_CaF2
				print 'MODO: %i' % modo

			elif flag_todo == 1:
				KBr_ensamble = KBr_en	### ENSAMBLE DE KBR CON VALORES DE SZA SIMILARES A CAF2

			KBr_sorted_idx = np.argsort(KBr_ensamble['epochTime'])
			KBr_ensamble = KBr_ensamble[KBr_sorted_idx]
			np.savetxt('CaF2ensamble.txt',CaF2_ensamble,fmt='%s %i %.04f %.04f %.04e %.04e %.04f %.04f %.04f %.04f')
			np.savetxt('KBrensamble%i.txt' % flag_todo,KBr_ensamble,fmt='%s %i %.04f %.04f %.04e %.04e %.04f %.04f %.04f %.04f')
			print 'ARCHIVOS CaF2ensamble.txt y KBrensamble%i.txt GENERADOS ;)' %flag_todo

		unq_arr_KBr,unique_idx,inverse_KBr = np.unique(KBr_ensamble['fname'],return_index=True,return_inverse=True)	
		print 'REPETIDOS (TOTALES: %i): %i' % (len(KBr_ensamble['fname']),len(KBr_ensamble['fname']) - len(unq_arr_KBr))
		if len(KBr_ensamble['fname']) - len(unq_arr_KBr) != 0:
			print 'ARCHIVOS REPETIDOS: '
			for element in inverse_KBr:
				if element not in unique_idx:
					print KBr_ensamble['fname'][element],KBr_ensamble['SZA'][element]
				
		KBr_f_dt = np.array([dt.datetime.fromtimestamp(uu) for uu in KBr_ensamble['epochTime']])
		KBr_f_d = np.array([uu.date() for uu in KBr_f_dt])
		CaF2_f_dt = np.array([dt.datetime.fromtimestamp(uu) for uu in CaF2_ensamble['epochTime']])
		CaF2_f_d = np.array([uu.date() for uu in CaF2_f_dt])

		unq_d_KBr = np.unique(KBr_f_d)
		for ele in unq_d_KBr:
			sel_KBr_dt = KBr_f_dt[KBr_f_d == ele]
			sel_KBr_dt = np.sort(sel_KBr_dt)
			#print 'KBr',ele,len(sel_KBr_dt),sel_KBr_dt[0].time(),sel_KBr_dt[-1].time(),'TIME DIF',(sel_KBr_dt[-1]-sel_KBr_dt[0]).total_seconds()/3600.0

		unq_d_CaF2 = np.unique(CaF2_f_d)
		for ele in unq_d_CaF2:
			sel_CaF2_dt = CaF2_f_dt[CaF2_f_d == ele]
			sel_CaF2_dt = np.sort(sel_CaF2_dt)
			#print 'CaF2',ele,len(sel_CaF2_dt),sel_CaF2_dt[0].time(),sel_CaF2_dt[-1].time(),'TIME DIF',(sel_CaF2_dt[-1]-sel_CaF2_dt[0]).total_seconds()/3600.0
	
		ens_hist_SZA_KBr, ens_be_SZA_KBr = np.histogram(KBr_ensamble['SZA'],be_SZA_KBr)
		ens_hist_SZA_norm_KBr = np.array([float(ele)/float(np.sum(ens_hist_SZA_KBr)) for ele in ens_hist_SZA_KBr])
		ens_bp_SZA_KBr = np.array([(ele + ens_be_SZA_KBr[ii+1])/2 for ii,ele in enumerate(ens_be_SZA_KBr[:-1])])

		mn_CO2_CaF2 = np.mean(CaF2_ensamble['totcolCO2'])
		std_CO2_CaF2 = np.std(CaF2_ensamble['totcolCO2'])
		mn_CO2_KBr = np.mean(KBr_ensamble['totcolCO2'])
		std_CO2_KBr = np.std(KBr_ensamble['totcolCO2'])
		mn_O2_CaF2 = np.mean(CaF2_ensamble['totcolO2'])
		std_O2_CaF2 = np.std(CaF2_ensamble['totcolO2'])
		mn_O2_KBr = np.mean(KBr_ensamble['totcolO2'])
		std_O2_KBr = np.std(KBr_ensamble['totcolO2'])
		mn_XCO2_CaF2 = np.mean(CaF2_ensamble['XCO2'])
		std_XCO2_CaF2 = np.std(CaF2_ensamble['XCO2'])
		mn_XCO2_KBr = np.mean(KBr_ensamble['XCO2'])
		std_XCO2_KBr = np.std(KBr_ensamble['XCO2'])

		hist_SZA_CaF2_all = np.array([ele + hist_SZA_CaF2_2015[ii] for ii,ele in enumerate(hist_SZA_CaF2_2014)])
		hist_SZA_norm_CaF2_all = np.array([float(ele)/float(np.sum(hist_SZA_CaF2_all)) for ele in hist_SZA_CaF2_all])

		plt.figure()
		#plt.suptitle('SZA ('+archivo_txt+')')
		#plt.plot(bp_SZA_CaF2_2014,hist_SZA_norm_CaF2_2014,'--b',label='CaF2 2014')
		#plt.plot(bp_SZA_CaF2_2015,hist_SZA_norm_CaF2_2015,'-.b',label='CaF2 2015')
		plt.plot(bp_SZA_CaF2_2015,hist_SZA_norm_CaF2_all,color='b',linestyle='-',label='CaF2 ENSAMBLE',linewidth=5)
		plt.plot(bp_SZA_KBr,hist_SZA_norm_KBr,color='r',linestyle=':',label='KBr ALL',linewidth=3)
		plt.plot(ens_bp_SZA_KBr,ens_hist_SZA_norm_KBr,color='r',linestyle='-.',label='KBr ENSAMBLE',linewidth=5	)

		plt.xlabel('SZA')
		plt.ylabel('Relative Frecuency')
		plt.legend(fontsize=12)
		#plt.savefig('SZA_DS.eps')

		#fecfec_KBr = np.array([dt.datetime.fromtimestamp(ele).replace(year=2001) for ele in KBr_ensamble['epochTime']])
		#fecfec_CaF2 = np.array([dt.datetime.fromtimestamp(ele).replace(year=2001) for ele in CaF2_ensamble['epochTime']])
		fecfec_KBr = np.array([dt.datetime.fromtimestamp(ele) for ele in KBr_ensamble['epochTime']])
		fecfec_CaF2 = np.array([dt.datetime.fromtimestamp(ele) for ele in CaF2_ensamble['epochTime']])


		plt.figure()
		#plt.suptitle('COLUMNS ('+archivo_txt+')')
		plt.subplot(311)
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
		#plt.xticks(rotation=45,fontsize=8)
		plt.xticks(fontsize=6)
		plt.plot(fecfec_KBr,KBr_ensamble['totcolCO2']/(10**21),'.r',label='KBr (%i)' % (len(KBr_ensamble['totcolCO2'])))
		plt.plot(fecfec_CaF2,CaF2_ensamble['totcolCO2']/(10**21),'.b',label='CaF2 (%i)' % (len(CaF2_ensamble['totcolCO2'])))
		plt.ylabel('$CO_2\/(10^{21}\/molec/cm^2)$')
		plt.legend(fontsize=8,ncol=2,loc=4)
		plt.subplot(312)
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
		#plt.xticks(rotation=45,fontsize=8)
		plt.xticks(fontsize=6)
		plt.plot(fecfec_KBr,KBr_ensamble['totcolO2']/(10**24),'.r',label='KBr')
		plt.plot(fecfec_CaF2,CaF2_ensamble['totcolO2']/(10**24),'.b',label='CaF2')
		plt.ylabel('$O_2\/(10^{24}\/molec/cm^2)$')
		#plt.legend()
		plt.subplot(313)
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
		#plt.xticks(rotation=45,fontsize=8)
		plt.xticks(fontsize=6)
		plt.plot(fecfec_KBr,KBr_ensamble['XCO2'],'.r',label='KBr')
		plt.plot(fecfec_CaF2,CaF2_ensamble['XCO2'],'.b',label='CaF2')
		#plt.xlabel('Epoch Time')
		plt.ylabel('$XCO_2\/(ppm)$')
		#plt.legend()
		if flag_todo == 1:
			pass
			#plt.savefig('COLUMNS_ENS.jpg')
		else:		
			pass
			#plt.savefig('COLUMNS_ENSII.jpg')

		return KBr_ensamble,CaF2_ensamble		

	#############################
	
	def HOUR_DTR(archivo_txt,f_col,KBr_idx,CaF2_idx):		### PARA UN SOLO DIA
	
		f_dt = np.array([dt.datetime.fromtimestamp(ele) for ele in f_col['epochTime']])
		f_hour = np.array([dt.datetime.fromtimestamp(ele).time().hour for ele in f_col['epochTime']])

		# HOURLY MEANS
		unq_hour = np.unique(f_hour)
		hourlymean_CO2 = np.array([np.mean(f_col['totcolCO2'][f_hour == ele]) for ele in unq_hour])
		hourlymean_O2 = np.array([np.mean(f_col['totcolO2'][f_hour == ele]) for ele in unq_hour])
		hourlymean_XCO2 = np.array([np.mean(f_col['XCO2'][f_hour == ele]) for ele in unq_hour])
		mean_hour = np.array([dt.datetime(2015,11,11,ele,30) for ele in unq_hour])

		#DETRENDING WITH HOURLY MEAN
		dt_dc = np.dtype([('epochTime',int),('dtr_CO2',float),('dtr_O2',float),('dtr_XCO2',float)])
		dtr_hr = np.empty([len(f_col['totcolCO2'])],dtype=dt_dc)
		dtr_hr['epochTime'] = f_col['epochTime']
		for ii,ele in enumerate(f_hour):
			dtr_hr['dtr_CO2'][ii] = f_col['totcolCO2'][ii] - hourlymean_CO2[unq_hour == ele]
			dtr_hr['dtr_O2'][ii] = f_col['totcolO2'][ii] - hourlymean_O2[unq_hour == ele]
			dtr_hr['dtr_XCO2'][ii] = f_col['XCO2'][ii] - hourlymean_XCO2[unq_hour == ele]

		plt.figure()
		#plt.suptitle('HOURLY MEANS ('+archivo_txt+')')
		plt.subplot(311)
		plt.plot(f_dt[CaF2_idx],f_col['totcolCO2'][CaF2_idx],'.b',label='CaF2')
		plt.plot(f_dt[KBr_idx],f_col['totcolCO2'][KBr_idx],'.r',label='KBr')
		#plt.plot(mean_hour,hourlymean_CO2,'--*g')
		plt.ylabel('CO2 (molecules/cm2)')
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
		plt.xticks(rotation=45,fontsize=8)
		#plt.legend(fontsize=8,loc=8,ncol=2)		
		plt.subplot(312)
		plt.plot(f_dt[CaF2_idx],f_col['totcolO2'][CaF2_idx],'.b',label='CaF2')
		plt.plot(f_dt[KBr_idx],f_col['totcolO2'][KBr_idx],'.r',label='KBr')
		#plt.plot(mean_hour,hourlymean_O2,'--*g')
		plt.ylabel('O2 (molecules/cm2)')
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
		plt.xticks(rotation=45,fontsize=8)
		plt.legend(fontsize=8,loc=8,ncol=2)		
		plt.subplot(313)
		plt.plot(f_dt[CaF2_idx],f_col['XCO2'][CaF2_idx],'.b',label='CaF2')
		plt.plot(f_dt[KBr_idx],f_col['XCO2'][KBr_idx],'.r',label='KBr')
		#plt.plot(mean_hour,hourlymean_XCO2,'--*g')
		plt.ylabel('XCO2 (ppm)')
		plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
		plt.xticks(rotation=45,fontsize=8)
		#plt.legend(fontsize=8,loc=8,ncol=2)
		#plt.savefig('COLUMNS_HR.eps')

		return dtr_hr

	##############################

	def DAY_DTR(archivo_txt,f_col,KBr_idx,CaF2_idx):			### NO SIRVE, REVISAR 05/02/2016

		f_date = np.array([dt.datetime.fromtimestamp(ele).date() for ele in f_col['epochTime']])

		# DAILY MEANS
		unq_date = np.unique(f_date)
		dailymean_CO2 = np.array([np.mean(f_col['totcolCO2'][f_date == ele]) for ele in unq_date])
		dailymean_O2 = np.array([np.mean(f_col['totcolO2'][f_date == ele]) for ele in unq_date])
		dailymean_XCO2 = np.array([np.mean(f_col['XCO2'][f_date == ele]) for ele in unq_date])

		# DETRENDING WITH DAILY MEAN
		dtr_CO2 = np.array([f_col['totcolCO2'][ii] - dailymean_CO2[unq_date == ele] for ii,ele in enumerate(f_date)])
		dtr_O2 = np.array([f_col['totcolO2'][ii] - dailymean_O2[unq_date == ele] for ii,ele in enumerate(f_date)])
		dtr_XCO2 = np.array([f_col['XCO2'][ii] - dailymean_XCO2[unq_date == ele] for ii,ele in enumerate(f_date)])

		f_fechasCaF2 = np.genfromtxt('fechasCaF2.txt',dtype=str)
		f_fechasCaF2 = np.array([dt.datetime.strptime(ele,'%y%m%d').date() for ele in f_fechasCaF2])
		#print f_fechasCaF2[0:10]

		CaF2_idx = np.in1d(f_date,f_fechasCaF2)
		print 'CaF2',len(f_col['totcolCO2'][CaF2_idx])
		CaF2_idx_dm = np.in1d(unq_date,f_fechasCaF2)		### SOLO PARA PROMEDIOS DIARIOS
		KBr_idx = np.invert(CaF2_idx)
		print 'KBr',len(f_col['totcolCO2'][KBr_idx])
		KBr_idx_dm = np.invert(CaF2_idx_dm)			### SOLO PARA PROMEDIOS DIARIOS

		#plt.figure()
		#plt.suptitle('DAILY MEANS ('+archivo_txt+')')
		#plt.subplot(311)
		#plt.plot(unq_date[CaF2_idx_dm],dailymean_CO2[CaF2_idx_dm],'.b')
		#plt.plot(unq_date[KBr_idx_dm],dailymean_CO2[KBr_idx_dm],'.r')
		#plt.ylabel('CO2 (molecules/cm2)')
		#plt.subplot(312)
		#plt.plot(unq_date[CaF2_idx_dm],dailymean_O2[CaF2_idx_dm],'.b')
		#plt.plot(unq_date[KBr_idx_dm],dailymean_O2[KBr_idx_dm],'.r')
		#plt.ylabel('O2 (molecules/cm2)')
		#plt.subplot(313)
		#plt.plot(unq_date[CaF2_idx_dm],dailymean_XCO2[CaF2_idx_dm],'.b')
		#plt.plot(unq_date[KBr_idx_dm],dailymean_XCO2[KBr_idx_dm],'.r')
		#plt.ylabel('XCO2 (ppm)')

	##############################
	
	def CONT_DTR (ensamble,bms):

		# DETRENDING OF CONTIGUOUS MEASUREMENTS
		f_date = np.array([dt.datetime.fromtimestamp(ele).date() for ele in ensamble['epochTime']])
		f_time = np.array([dt.datetime.fromtimestamp(ele).time() for ele in ensamble['epochTime']])
		unq_date = np.unique(f_date)
		suma_CO2 = []
		suma_O2 = []
		suma_XCO2 = []
		hrs = []
		dtr_cont_CO2 = []
		dtr_cont_O2 = []
		dtr_cont_XCO2 = []
		std_cont_CO2 = []
		std_cont_O2 = []
		std_cont_XCO2 = []
		for unq_d in unq_date:
			fh = f_time[f_date == unq_d]
			ch_CO2 = ensamble['totcolCO2'][f_date == unq_d]
			ch_O2 = ensamble['totcolO2'][f_date == unq_d]
			ch_XCO2 = ensamble['XCO2'][f_date == unq_d]
			idx_fh = fh.argsort()
			ch_CO2 = ch_CO2[idx_fh]
			ch_O2 = ch_O2[idx_fh]
			ch_XCO2 = ch_XCO2[idx_fh]
			fh = fh[idx_fh]
			#print unq_d
			#print fh
			datehour = np.array([dt.datetime.combine(unq_d,ele) for ele in fh])
			cont_6 = 0
			for ii,hora in enumerate(datehour):
				if ii != 0:
					if (((hora - datehour[ii-1]).total_seconds() > 180.0) | (cont_6 > 5)):
						mn_CO2 = np.mean(suma_CO2)
						mn_O2 = np.mean(suma_O2)
						mn_XCO2 = np.mean(suma_XCO2)
						suma_CO2 = np.array([100.0*(ele - mn_CO2)/mn_CO2 for ele in suma_CO2])
						suma_O2 = np.array([100.0*(ele - mn_O2)/mn_O2 for ele in suma_O2])
						suma_XCO2 = np.array([100.0*(ele - mn_XCO2)/mn_XCO2 for ele in suma_XCO2])
						dtr_cont_CO2.extend(suma_CO2)
						dtr_cont_O2.extend(suma_O2)
						dtr_cont_XCO2.extend(suma_XCO2)
						if len(suma_CO2) > 1:
							std_cont_CO2.append(np.std(suma_CO2,ddof = 1))
							std_cont_O2.append(np.std(suma_O2,ddof = 1))
							std_cont_XCO2.append(np.std(suma_XCO2,ddof = 1))
						#else:
						#	print hora,hrs
						print >> open('DESV_ENS.txt','a'), 'FECHA (%s): %s' % (bms,unq_d)
						print >> open('DESV_ENS.txt','a'), [str(ele.time()) for ele in hrs]
						print >> open('DESV_ENS.txt','a'), 'STD CO2: %.4f, STD O2: %.4f, STD XCO2: %.4f' % (np.std(suma_CO2),np.std(suma_O2),np.std(suma_XCO2))
						suma_CO2 = []
						suma_O2 = []
						suma_XCO2 = []
						hrs = []
						cont_6 = 0
				suma_CO2.append(ch_CO2[ii])
				suma_O2.append(ch_O2[ii])
				suma_XCO2.append(ch_XCO2[ii])
				hrs.append(hora)
				cont_6 = cont_6 + 1
			mn_CO2 = np.mean(suma_CO2)
			mn_O2 = np.mean(suma_O2)
			mn_XCO2 = np.mean(suma_XCO2)
			suma_CO2 = np.array([(ele - mn_CO2)/mn_CO2 for ele in suma_CO2])
			suma_O2 = np.array([(ele - mn_O2)/mn_O2 for ele in suma_O2])
			suma_XCO2 = np.array([(ele - mn_XCO2)/mn_XCO2 for ele in suma_XCO2])
			dtr_cont_CO2.extend(suma_CO2)
			dtr_cont_O2.extend(suma_O2)
			dtr_cont_XCO2.extend(suma_XCO2)
			if len(suma_CO2) > 1:
				std_cont_CO2.append(np.std(suma_CO2,ddof = 1))
				std_cont_O2.append(np.std(suma_O2,ddof = 1))
				std_cont_XCO2.append(np.std(suma_XCO2,ddof = 1))
			else:
				print bms,hora,hrs
			#print unq_d
			#print 'HORAS',hrs
			#print suma_CO2
			suma_CO2 = []
			suma_O2 = []
			suma_XCO2 = []
			hrs = []
			cont_6 = 0

		print 'NUMBER OF GROUPS OF CONTIGUOUS MEAS (%s,TOTAL: %i): %i' % (bms,len(dtr_cont_CO2),len(std_cont_CO2))
		dtr_cont_CO2 = np.array(dtr_cont_CO2)
		dtr_cont_O2 = np.array(dtr_cont_O2)
		dtr_cont_XCO2 = np.array(dtr_cont_XCO2)
		std_cont_CO2 = np.array(std_cont_CO2)
		std_cont_O2 = np.array(std_cont_O2)
		std_cont_XCO2 = np.array(std_cont_XCO2)
		
		dt_dc = np.dtype([('epochTime',int),('dtr_CO2',float),('dtr_O2',float),('dtr_XCO2',float)])

		dtr_cont = np.empty([len(dtr_cont_CO2)],dtype=dt_dc)
		dtr_cont['epochTime'] = ensamble['epochTime']
		dtr_cont['dtr_CO2'] = dtr_cont_CO2
		dtr_cont['dtr_O2'] = dtr_cont_O2
		dtr_cont['dtr_XCO2'] = dtr_cont_XCO2

		std_cont = np.empty([len(std_cont_CO2)],dtype=[('std_CO2',float),('std_O2',float),('std_XCO2',float)])
		std_cont['std_CO2'] = std_cont_CO2
		std_cont['std_O2'] = std_cont_O2
		std_cont['std_XCO2'] = std_cont_XCO2

		return dtr_cont,std_cont

	#print 'CaF2 CO2	',len(f_col['totcolCO2'][CaF2_idx]),np.mean(f_col['totcolCO2'][CaF2_idx]),np.std(f_col['totcolCO2'][CaF2_idx])
	#print 'CaF2 O2		',len(f_col['totcolO2'][CaF2_idx]),np.mean(f_col['totcolO2'][CaF2_idx]),np.std(f_col['totcolO2'][CaF2_idx])
	#print 'CaF2 XCO2	',len(f_col['XCO2'][CaF2_idx]),np.mean(f_col['XCO2'][CaF2_idx]),np.std(f_col['XCO2'][CaF2_idx])
	#print 'KBr CO2		',len(f_col['totcolCO2'][KBr_idx]),np.mean(f_col['totcolCO2'][KBr_idx]),np.std(f_col['totcolCO2'][KBr_idx])
	#print 'KBr O2		',len(f_col['totcolO2'][KBr_idx]),np.mean(f_col['totcolO2'][KBr_idx]),np.std(f_col['totcolO2'][KBr_idx])
	#print 'KBr XCO2	',len(f_col['XCO2'][KBr_idx]),np.mean(f_col['XCO2'][KBr_idx]),np.std(f_col['XCO2'][KBr_idx])

	#print 'DTR CaF2 CO2	',len(dtr_CO2[CaF2_idx]),np.mean(dtr_CO2[CaF2_idx]),np.std(dtr_CO2[CaF2_idx])
	#print 'DTR CaF2 O2	',len(dtr_O2[CaF2_idx]),np.mean(dtr_O2[CaF2_idx]),np.std(dtr_O2[CaF2_idx])
	#print 'DTR CaF2 XCO2	',len(dtr_XCO2[CaF2_idx]),np.mean(dtr_XCO2[CaF2_idx]),np.std(dtr_XCO2[CaF2_idx])
	#print 'DTR KBr CO2	',len(dtr_CO2[KBr_idx]),np.mean(dtr_CO2[KBr_idx]),np.std(dtr_CO2[KBr_idx])
	#print 'DTR KBr O2	',len(dtr_O2[KBr_idx]),np.mean(dtr_O2[KBr_idx]),np.std(dtr_O2[KBr_idx])
	#print 'DTR KBr XCO2	',len(dtr_XCO2[KBr_idx]),np.mean(dtr_XCO2[KBr_idx]),np.std(dtr_XCO2[KBr_idx])

	#print 'DTR CONTIGUOUS CaF2 CO2	',len(dtr_cont_CO2[CaF2_idx]),np.mean(dtr_cont_CO2[CaF2_idx]),np.std(dtr_cont_CO2[CaF2_idx])**2
	#print 'DTR CONTIGUOUS CaF2 O2	',len(dtr_cont_O2[CaF2_idx]),np.mean(dtr_cont_O2[CaF2_idx]),np.std(dtr_cont_O2[CaF2_idx])**2
	#print 'DTR CONTIGUOUS CaF2 XCO2	',len(dtr_cont_XCO2[CaF2_idx]),np.mean(dtr_cont_XCO2[CaF2_idx]),np.std(dtr_cont_XCO2[CaF2_idx])**2
	#print 'DTR CONTIGUOUS KBr CO2	',len(dtr_cont_CO2[KBr_idx]),np.mean(dtr_cont_CO2[KBr_idx]),np.std(dtr_cont_CO2[KBr_idx])**2
	#print 'DTR CONTIGUOUS KBr O2	',len(dtr_cont_O2[KBr_idx]),np.mean(dtr_cont_O2[KBr_idx]),np.std(dtr_cont_O2[KBr_idx])**2
	#print 'DTR CONTIGUOUS KBr XCO2	',len(dtr_cont_XCO2[KBr_idx]),np.mean(dtr_cont_XCO2[KBr_idx]),np.std(dtr_cont_XCO2[KBr_idx])**2

	### AQUI EMPIEZA COLUMNS !!!
	
	print "------------------------------"
	print 'ARCHIVO USADO:',archivo_txt

	dt_txt = np.dtype([('fname',np.str_,30),('epochTime',int),('SZA',float), \
		('XCO2',float),('totcolCO2',float),('totcolO2',float),('LOSerrCO2',float),('NOISEerrCO2',float), \
		('LOSerrO2',float),('NOISEerrO2',float)])

	f_col = np.genfromtxt(archivo_txt,dtype = dt_txt,delimiter = " ")
		
	### FILTER
	print 'ORIGINAL: %i' % len(f_col['XCO2'])
	fac_fil_CO2 = 2.0
	up_fil_CO2 = np.mean(f_col['XCO2']) + fac_fil_CO2*np.std(f_col['XCO2'])
	dw_fil_CO2 = np.mean(f_col['XCO2']) - fac_fil_CO2*np.std(f_col['XCO2'])
	con_CO2 = ((f_col['XCO2'] > dw_fil_CO2) & (f_col['XCO2'] < up_fil_CO2))
	f_col = f_col[con_CO2]
	print 'AFTER FILTER: %i' % len(f_col['XCO2'])
	#########

	f_col = f_col[f_col['epochTime'].argsort()]	### SORTING BASED ON EPOCH TIME
	f_date = np.array([dt.datetime.fromtimestamp(ele).date() for ele in f_col['epochTime']])
	print 'TOTAL DAYS: %i' % len(np.unique(f_date))

	f_fechasCaF2 = np.genfromtxt('fechasCaF2.txt',dtype=str)
	f_fechasCaF2 = np.array([dt.datetime.strptime(ele,'%y%m%d').date() for ele in f_fechasCaF2])
	#print f_fechasCaF2[0:10]

	CaF2_idx = np.in1d(f_date,f_fechasCaF2)
	#CaF2_idx_dm = np.in1d(unq_date,f_fechasCaF2)		### SOLO PARA PROMEDIOS DIARIOS
	KBr_idx = np.invert(CaF2_idx)
	#KBr_idx_dm = np.invert(CaF2_idx_dm)			### SOLO PARA PROMEDIOS DIARIOS

	#HISTOS_COL(archivo_txt,f_col,KBr_idx,CaF2_idx,numbins=50)
	KBr_ens,CaF2_ens = ENSAMBLE(archivo_txt,f_col,KBr_idx,CaF2_idx,flag_todo=0,modo=0)	### MODO: 1: FECHAS SIMILARES EN DIFERENTES ANOS, 0: FECHAS MAS RECIENTES
	### flag_todo: 0 PARA USAR UN ENSAMBLE IGUAL AL DE CaF2, 1 PARA USAR TODOS LOS DATOS DE KBr
	dtr_cont_KBr,std_cont_KBr = CONT_DTR(KBr_ens,'KBr')
	dtr_cont_CaF2,std_cont_CaF2 = CONT_DTR(CaF2_ens,'CaF2')
	
	dt_KBr = np.array([dt.datetime.fromtimestamp(ele) for ele in dtr_cont_KBr['epochTime']])
	d_KBr = np.array([dt.datetime.fromtimestamp(ele).date() for ele in dtr_cont_KBr['epochTime']])
	t_KBr = np.array([dt.datetime.combine(dt.date(2000,1,1),dt.datetime.fromtimestamp(ele).time()) for ele in dtr_cont_KBr['epochTime']])
	dt_CaF2 = np.array([dt.datetime.fromtimestamp(ele) for ele in dtr_cont_CaF2['epochTime']])
	d_CaF2 = np.array([dt.datetime.fromtimestamp(ele).date() for ele in dtr_cont_CaF2['epochTime']])
	t_CaF2 = np.array([dt.datetime.combine(dt.date(2000,1,1),dt.datetime.fromtimestamp(ele).time()) for ele in dtr_cont_CaF2['epochTime']])

	print 'DAYS KBr: %i, CaF2: %i' % (len(np.unique(d_KBr)),len(np.unique(d_CaF2)))

	mn_CO2_CaF2 = np.mean(CaF2_ens['totcolCO2'])
	std_CO2_CaF2 = np.std(CaF2_ens['totcolCO2'])
	mn_CO2_KBr = np.mean(KBr_ens['totcolCO2'])
	std_CO2_KBr = np.std(KBr_ens['totcolCO2'])
	mn_O2_CaF2 = np.mean(CaF2_ens['totcolO2'])
	std_O2_CaF2 = np.std(CaF2_ens['totcolO2'])
	mn_O2_KBr = np.mean(KBr_ens['totcolO2'])
	std_O2_KBr = np.std(KBr_ens['totcolO2'])
	mn_XCO2_CaF2 = np.mean(CaF2_ens['XCO2'])
	std_XCO2_CaF2 = np.std(CaF2_ens['XCO2'])
	mn_XCO2_KBr = np.mean(KBr_ens['XCO2'])
	std_XCO2_KBr = np.std(KBr_ens['XCO2'])

	print "------------------------------"
	print "ENSAMBLES KBr:",len(KBr_ens),'ELEMENTOS','CaF2:',len(CaF2_ens),'ELEMENTOS'
	'''
	print "CO2 (KBr)","%.6f %%" % ((100.0*std_CO2_KBr)/mn_CO2_KBr),"	O2 (KBr)","%.6f %%" % ((100.0*std_O2_KBr)/mn_O2_KBr),\
	"	XCO2 (KBr)","%.6f %%" % ((100.0*std_XCO2_KBr)/mn_XCO2_KBr)
	#print "CO2 (KBr)","%.6f %%" % ((100.0*np.std(f_col['totcolCO2'][KBr_idx]))/np.mean(f_col['totcolCO2'][KBr_idx])), \
	"O2 (KBr)","%.6f %%" % ((100.0*np.std(f_col['totcolO2'][KBr_idx]))/np.mean(f_col['totcolO2'][KBr_idx])), \
	"XCO2 (KBr)","%.6f %%" % ((100.0*np.std(f_col['XCO2'][KBr_idx]))/np.mean(f_col['XCO2'][KBr_idx]))	### TODOS LOS DATOS DE KBR
	
	print "CO2 (CaF2)","%.6f %%" % ((100.0*std_CO2_CaF2)/mn_CO2_CaF2),"	O2 (CaF2)","%.6f %%" % ((100.0*std_O2_CaF2)/mn_O2_CaF2),\
	"	XCO2 (CaF2)","%.6f %% \n" % ((100.0*std_XCO2_CaF2)/mn_XCO2_CaF2)

	#print 'SZA KBr 20',KBr_ens[(KBr_ens['SZA'] > 19.0) & (KBr_ens['SZA'] < 21.0)][0],'70',KBr_ens[(KBr_ens['SZA'] > 69.0) & (KBr_ens['SZA'] < 71.0)][0]
	#print 'SZA CaF2 20',CaF2_ens[(CaF2_ens['SZA'] > 19.0) & (CaF2_ens['SZA'] < 21.0)][0],'70',CaF2_ens[(CaF2_ens['SZA'] > 69.0) & (CaF2_ens['SZA'] < 71.0)][0]
	print "PRF LOS CO2 (KBr) %.6f PRF LOS O2 (KBr) %.6f PRF LOS TOTAL (KBr) %.6f" % (np.mean(KBr_ens['LOSerrCO2'])/0.9,np.mean(KBr_ens['LOSerrO2'])/0.9,np.sqrt((np.mean(KBr_ens['LOSerrCO2'])/0.9)**2 + (np.mean(KBr_ens['LOSerrO2'])/0.9)**2))
	print "PRF LOS CO2 (CaF2) %.6f PRF LOS O2 (CaF2) %.6f PRF LOS TOTAL (CaF2) %.6f" % (np.mean(CaF2_ens['LOSerrCO2'])/0.9,np.mean(CaF2_ens['LOSerrO2'])/0.9,np.sqrt((np.mean(CaF2_ens['LOSerrCO2'])/0.9)**2 + (np.mean(CaF2_ens['LOSerrO2'])/0.9)**2))

	print "PRF NOISE CO2 (KBr) %.6f	PRF NOISE O2 (KBr) %.6f PRF NOISE XCO2 (KBr) %.6f" % (np.mean(KBr_ens['NOISEerrCO2']),np.mean(KBr_ens['NOISEerrO2']),np.sqrt((np.mean(KBr_ens['NOISEerrCO2']))**2 + (np.mean(KBr_ens['NOISEerrO2'])**2)))
	print "PRF NOISE CO2 (CaF2) %.6f PRF NOISE O2 (CaF2) %.6f PRF NOISE XCO2 (CaF2) %.6f" % (np.mean(CaF2_ens['NOISEerrCO2']),np.mean(CaF2_ens['NOISEerrO2']),np.sqrt((np.mean(CaF2_ens['NOISEerrCO2'])**2) + (np.mean(CaF2_ens['NOISEerrO2'])**2)))
	print "------------------------------"
	'''
	sigma_O2_KBr = np.std(dtr_cont_KBr['dtr_O2'])#/np.mean(f_col['totcolO2'][KBr_idx])
	sigma_O2_CaF2 = np.std(dtr_cont_CaF2['dtr_O2'])#/np.mean(f_col['totcolO2'][CaF2_idx])
	sigma_CO2_KBr = np.std(dtr_cont_KBr['dtr_CO2'])#/np.mean(f_col['totcolCO2'][KBr_idx])
	sigma_CO2_CaF2 = np.std(dtr_cont_CaF2['dtr_CO2'])#/np.mean(f_col['totcolCO2'][CaF2_idx])
	sigma_XCO2_KBr = np.std(dtr_cont_KBr['dtr_XCO2'])#/np.mean(f_col['XCO2'][KBr_idx])
	sigma_XCO2_CaF2 = np.std(dtr_cont_CaF2['dtr_XCO2'])#/np.mean(f_col['XCO2'][CaF2_idx])

	sigma_O2_KBrII = np.std(dtr_cont_KBr['dtr_O2'][dtr_cont_KBr['dtr_O2'] != 0])#/np.mean(f_col['totcolO2'][KBr_idx])
	sigma_O2_CaF2II = np.std(dtr_cont_CaF2['dtr_O2'][dtr_cont_CaF2['dtr_O2'] != 0])#/np.mean(f_col['totcolO2'][CaF2_idx])
	sigma_CO2_KBrII = np.std(dtr_cont_KBr['dtr_CO2'][dtr_cont_KBr['dtr_CO2'] != 0])#/np.mean(f_col['totcolCO2'][KBr_idx])
	sigma_CO2_CaF2II = np.std(dtr_cont_CaF2['dtr_CO2'][dtr_cont_CaF2['dtr_CO2'] != 0])#/np.mean(f_col['totcolCO2'][CaF2_idx])
	sigma_XCO2_KBrII = np.std(dtr_cont_KBr['dtr_XCO2'][dtr_cont_KBr['dtr_XCO2'] != 0])#/np.mean(f_col['XCO2'][KBr_idx])
	sigma_XCO2_CaF2II = np.std(dtr_cont_CaF2['dtr_XCO2'][dtr_cont_CaF2['dtr_XCO2'] != 0])#/np.mean(f_col['XCO2'][CaF2_idx])

	sigma_O2_KBrIII = np.mean(std_cont_KBr['std_O2'])#/np.mean(f_col['totcolO2'][KBr_idx])
	sigma_O2_CaF2III = np.mean(std_cont_CaF2['std_O2'])#/np.mean(f_col['totcolO2'][CaF2_idx])
	sigma_CO2_KBrIII = np.mean(std_cont_KBr['std_CO2'])#/np.mean(f_col['totcolCO2'][KBr_idx])
	sigma_CO2_CaF2III = np.mean(std_cont_CaF2['std_CO2'])#/np.mean(f_col['totcolCO2'][CaF2_idx])
	sigma_XCO2_KBrIII = np.mean(std_cont_KBr['std_XCO2'])#/np.mean(f_col['XCO2'][KBr_idx])
	sigma_XCO2_CaF2III = np.mean(std_cont_CaF2['std_XCO2'])#/np.mean(f_col['XCO2'][CaF2_idx])

	print "------------------------------"
	print "CONTIGUOUS DTR"
	print "KBr	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%" % (sigma_CO2_KBr,sigma_O2_KBr,sigma_XCO2_KBr)
	print "CaF2	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%" % (sigma_CO2_CaF2,sigma_O2_CaF2,sigma_XCO2_CaF2)

	#print 'ZEROS IN dtr_cont_KBr:',len(dtr_cont_KBr[dtr_cont_KBr['dtr_CO2'] == 0])
	#print 'ZEROS IN dtr_cont_CaF2:',len(dtr_cont_CaF2[dtr_cont_CaF2['dtr_CO2'] == 0]),'\n'
	print 'NO ZEROES'
	print "KBr	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%" % (sigma_CO2_KBrII,sigma_O2_KBrII,sigma_XCO2_KBrII)
	print "CaF2	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%" % (sigma_CO2_CaF2II,sigma_O2_CaF2II,sigma_XCO2_CaF2II)

	'''
	K = np.ones(shape=(5,5))
	K[0,2] = K[0,3] = K[0,4] = 0#K[0,5] = 0
	K[1,1] = K[1,3] = K[1,4] = 0
	K[2,1] = K[2,2] = K[2,4] = 0#K[2,5] = 0
	K[3,1] = K[3,2] = K[3,3] = 0
	#K[4,0] = K[4,2] = K[4,4] = 0#K[4,5] = 0
	K[4,0] = K[4,1] = K[4,3] = 0#K[4,5] = 0
	'''
	K = np.ones(shape=(3,3))
	K[0,2] = K[1,1] = K[2,0] = 0

	K = K.T
	Kinv = np.linalg.inv(K)

	vect_sig_KBr = [sigma_CO2_KBrII**2,sigma_O2_KBrII**2,sigma_XCO2_KBrII**2]#,sigma_XCO2_CaF2**2]
	vect_sig_CaF2 = [sigma_CO2_CaF2II**2,sigma_O2_CaF2II**2,sigma_XCO2_CaF2II**2]
	vect_err_KBr = np.dot(vect_sig_KBr,Kinv)
	vect_err_CaF2 = np.dot(vect_sig_CaF2,Kinv)
	vect_err_root_KBr = [np.sqrt(ele) for ele in vect_err_KBr]
	vect_err_root_CaF2 = [np.sqrt(ele) for ele in vect_err_CaF2]

	vect_sig_KBrIII = [sigma_CO2_KBrIII**2,sigma_O2_KBrIII**2,sigma_XCO2_KBrIII**2]#,sigma_XCO2_CaF2**2]
	vect_sig_CaF2III = [sigma_CO2_CaF2III**2,sigma_O2_CaF2III**2,sigma_XCO2_CaF2III**2]
	vect_err_KBrIII = np.dot(vect_sig_KBrIII,Kinv)
	vect_err_CaF2III = np.dot(vect_sig_CaF2III,Kinv)
	vect_err_root_KBrIII = [np.sqrt(ele) for ele in vect_err_KBrIII]
	vect_err_root_CaF2III = [np.sqrt(ele) for ele in vect_err_CaF2III]
		

	print "------------------------------"
	print "KBr	CORRELATED: %.6F %%	NOISE CO2: %.6F %%	NOISE O2 %.6F %%	TOTAL XCO2 %.6f %%" % (vect_err_root_KBr[0],vect_err_root_KBr[1],vect_err_root_KBr[2],np.sqrt(vect_err_root_KBr[1]**2 + vect_err_root_KBr[2]**2))
	print "CaF2	CORRELATED: %.6F %%	NOISE CO2: %.6F %%	NOISE O2 %.6F %%	TOTAL XCO2 %.6f %%" % (vect_err_root_CaF2[0],vect_err_root_CaF2[1],vect_err_root_CaF2[2],np.sqrt(vect_err_root_CaF2[1]**2 + vect_err_root_CaF2[2]**2))
	print "------------------------------"
	'''
	print 'DIFF'
	print 'CORR: %.6f%%	REP CO2: %.6f%%	REP O2: %.6f%%	NOISE CO2: %.6f%%	NOISE O2: %.6f%%	TOTAL XCO2: %.6f%%' % \
	(100*(1-(vect_err_root_CaF2[0]/vect_err_root_KBr[0])),100*(1-(sigma_CO2_CaF2II/sigma_CO2_KBrII)),100*(1-(sigma_O2_CaF2II/sigma_O2_KBrII)),100*(1-(vect_err_root_CaF2[1]/vect_err_root_KBr[1])),100*(1-(vect_err_root_CaF2[2]/vect_err_root_KBr[2])),100*(1-(sigma_XCO2_CaF2II/sigma_XCO2_KBrII)) )
	print "------------------------------"
	'''
	print 'CHECK'
	print 'KBr	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%' % (np.sqrt(vect_err_root_KBr[0]**2+vect_err_root_KBr[1]**2),np.sqrt(vect_err_root_KBr[0]**2+vect_err_root_KBr[2]**2),np.sqrt(vect_err_root_KBr[1]**2+vect_err_root_KBr[2]**2))
	print 'CaF2	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%' % (np.sqrt(vect_err_root_CaF2[0]**2+vect_err_root_CaF2[1]**2),np.sqrt(vect_err_root_CaF2[0]**2+vect_err_root_CaF2[2]**2),np.sqrt(vect_err_root_CaF2[1]**2+vect_err_root_CaF2[2]**2))
	print "------------------------------"
	print "------------------------------"
	print 'OTRA FORMA'
	print "KBr	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%" % (sigma_CO2_KBrIII,sigma_O2_KBrIII,sigma_XCO2_KBrIII)
	print "CaF2	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%" % (sigma_CO2_CaF2III,sigma_O2_CaF2III,sigma_XCO2_CaF2III)
	print "------------------------------"
	print "KBr	CORRELATED: %.6F %%	NOISE CO2: %.6F %%	NOISE O2 %.6F %%	TOTAL XCO2 %.6f %%" % (vect_err_root_KBrIII[0],vect_err_root_KBrIII[1],vect_err_root_KBrIII[2],np.sqrt(vect_err_root_KBrIII[1]**2 + vect_err_root_KBrIII[2]**2))
	print "CaF2	CORRELATED: %.6F %%	NOISE CO2: %.6F %%	NOISE O2 %.6F %%	TOTAL XCO2 %.6f %%" % (vect_err_root_CaF2III[0],vect_err_root_CaF2III[1],vect_err_root_CaF2III[2],np.sqrt(vect_err_root_CaF2III[1]**2 + vect_err_root_CaF2III[2]**2))
	print "------------------------------"
	'''
	print 'DIFF'
	print 'CORR: %.6f%%	REP CO2: %.6f%%	REP O2: %.6f%%	NOISE CO2: %.6f%%	NOISE O2: %.6f%%	TOTAL XCO2: %.6f%%' % \
	(100*(1-(vect_err_root_CaF2III[0]/vect_err_root_KBrIII[0])),100*(1-(sigma_CO2_CaF2III/sigma_CO2_KBrIII)),100*(1-(sigma_O2_CaF2III/sigma_O2_KBrIII)),100*(1-(vect_err_root_CaF2III[1]/vect_err_root_KBrIII[1])),100*(1-(vect_err_root_CaF2III[2]/vect_err_root_KBrIII[2])),100*(1-(sigma_XCO2_CaF2III/sigma_XCO2_KBrIII)) )
	print "------------------------------"
	'''
	print 'CHECK'
	print 'KBr	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%' % (np.sqrt(vect_err_root_KBrIII[0]**2+vect_err_root_KBrIII[1]**2),np.sqrt(vect_err_root_KBrIII[0]**2+vect_err_root_KBrIII[2]**2),np.sqrt(vect_err_root_KBrIII[1]**2+vect_err_root_KBrIII[2]**2))
	print 'CaF2	CO2: %.6f %%	O2: %.6f %%	XCO2: %.6f %%' % (np.sqrt(vect_err_root_CaF2III[0]**2+vect_err_root_CaF2III[1]**2),np.sqrt(vect_err_root_CaF2III[0]**2+vect_err_root_CaF2III[2]**2),np.sqrt(vect_err_root_CaF2III[1]**2+vect_err_root_CaF2III[2]**2))

	n_bins_hist = 35	
	hist_KBr_CO2,be_KBr_CO2 = np.histogram(std_cont_KBr['std_CO2'],bins=n_bins_hist)
	hist_KBr_CO2 = [float(ele)/np.sum(hist_KBr_CO2) for ele in hist_KBr_CO2]
	bp_KBr_CO2 = [(ele + be_KBr_CO2[ee])/2 for ee,ele in enumerate(be_KBr_CO2[:-1])]
	cum_func_KBr_CO2 = [hist_KBr_CO2[0]]
	for uu,ulu in enumerate(hist_KBr_CO2[1:]):
		cum_func_KBr_CO2.append(ulu + cum_func_KBr_CO2[uu])
	hist_KBr_O2,be_KBr_O2 = np.histogram(std_cont_KBr['std_O2'],bins=n_bins_hist)
	hist_KBr_O2 = [float(ele)/np.sum(hist_KBr_O2) for ele in hist_KBr_O2]
	bp_KBr_O2 = [(ele + be_KBr_O2[ee])/2 for ee,ele in enumerate(be_KBr_O2[:-1])]
	cum_func_KBr_O2 = [hist_KBr_O2[0]]
	for uu,ulu in enumerate(hist_KBr_O2[1:]):
		cum_func_KBr_O2.append(ulu + cum_func_KBr_O2[uu])
	hist_KBr_XCO2,be_KBr_XCO2 = np.histogram(std_cont_KBr['std_XCO2'],bins=n_bins_hist)
	hist_KBr_XCO2 = [float(ele)/np.sum(hist_KBr_XCO2) for ele in hist_KBr_XCO2]
	bp_KBr_XCO2 = [(ele + be_KBr_XCO2[ee])/2 for ee,ele in enumerate(be_KBr_XCO2[:-1])]
	cum_func_KBr_XCO2 = [hist_KBr_XCO2[0]]
	for uu,ulu in enumerate(hist_KBr_XCO2[1:]):
		cum_func_KBr_XCO2.append(ulu + cum_func_KBr_XCO2[uu])
	hist_CaF2_CO2,be_CaF2_CO2 = np.histogram(std_cont_CaF2['std_CO2'],bins=n_bins_hist)
	hist_CaF2_CO2 = [float(ele)/np.sum(hist_CaF2_CO2) for ele in hist_CaF2_CO2]
	bp_CaF2_CO2 = [(ele + be_CaF2_CO2[ee])/2 for ee,ele in enumerate(be_CaF2_CO2[:-1])]
	cum_func_CaF2_CO2 = [hist_CaF2_CO2[0]]
	for uu,ulu in enumerate(hist_CaF2_CO2[1:]):
		cum_func_CaF2_CO2.append(ulu + cum_func_CaF2_CO2[uu])
	hist_CaF2_O2,be_CaF2_O2 = np.histogram(std_cont_CaF2['std_O2'],bins=n_bins_hist)
	hist_CaF2_O2 = [float(ele)/np.sum(hist_CaF2_O2) for ele in hist_CaF2_O2]
	bp_CaF2_O2 = [(ele + be_CaF2_O2[ee])/2 for ee,ele in enumerate(be_CaF2_O2[:-1])]
	cum_func_CaF2_O2 = [hist_CaF2_O2[0]]
	for uu,ulu in enumerate(hist_CaF2_O2[1:]):
		cum_func_CaF2_O2.append(ulu + cum_func_CaF2_O2[uu])
	hist_CaF2_XCO2,be_CaF2_XCO2 = np.histogram(std_cont_CaF2['std_XCO2'],bins=n_bins_hist)
	hist_CaF2_XCO2 = [float(ele)/np.sum(hist_CaF2_XCO2) for ele in hist_CaF2_XCO2]
	bp_CaF2_XCO2 = [(ele + be_CaF2_XCO2[ee])/2 for ee,ele in enumerate(be_CaF2_XCO2[:-1])]
	cum_func_CaF2_XCO2 = [hist_CaF2_XCO2[0]]
	for ee,ele in enumerate(hist_CaF2_XCO2[1:]):
		cum_func_CaF2_XCO2.append(ele + cum_func_CaF2_XCO2[ee])

	plt.figure()
	plt.subplot(211)
	plt.title('CO2')
	plt.plot(bp_KBr_CO2,hist_KBr_CO2,color='r',label='KBr',linestyle='-')
	plt.plot(bp_CaF2_CO2,hist_CaF2_CO2,color='b',label='CaF2',linestyle='-')
	#plt.tick_params(axis='x',which='both', labelbottom='off')
	#plt.legend()
	plt.ylabel('Relative Frecuency')
	plt.subplot(212)
	plt.title('O2')
	plt.plot(bp_KBr_O2,hist_KBr_O2,color='r',label='KBr',linestyle='-')
	plt.plot(bp_CaF2_O2,hist_CaF2_O2,color='b',label='CaF2',linestyle='-')
	#plt.legend()
	plt.ylabel('Relative Frecuency')
	plt.xlabel('Standard Deviation')
	#plt.savefig('HIST_COLS.jpg')

	plt.figure()
	plt.suptitle('XCO2')
	plt.plot([0.25,0.25],[0.0,max(hist_KBr_XCO2)],linestyle='--',color='k')
	plt.plot([sigma_XCO2_KBrII,sigma_XCO2_KBrII],[0.0,max(hist_KBr_XCO2)],linestyle='--',color='r')
	plt.plot([sigma_XCO2_CaF2II,sigma_XCO2_CaF2II],[0.0,max(hist_KBr_XCO2)],linestyle='--',color='b')
	plt.plot(bp_KBr_XCO2,hist_KBr_XCO2,color='r',label='KBr',linestyle='-')
	plt.plot(bp_CaF2_XCO2,hist_CaF2_XCO2,color='b',label='CaF2',linestyle='-')
	plt.legend()
	plt.xlabel('Standard Deviation')
	plt.ylabel('Relative Frecuency')
	#plt.savefig('HIST_XCO2.jpg')

	conf = 0.98
	plt.figure()
	#plt.suptitle('$XCO_2$')
	#plt.plot([0.0,max(bp_KBr_XCO2)],[conf,conf],linestyle=':',color='k')
	plt.plot([0.25,0.25],[0.0,1.0],linestyle='--',color='k')
	plt.text(0.30,0.6,'TCCON\nPrecision',color='k',horizontalalignment='center')
	plt.plot([sigma_XCO2_KBrIII,sigma_XCO2_KBrIII],[0.0,1.0],linestyle='--',color='r')
	plt.text(0.14,0.6,'Mean\n'+r'$XCO_2$'+'\nPrecision',color='r',horizontalalignment='center')
	plt.plot([sigma_XCO2_CaF2III,sigma_XCO2_CaF2III],[0.0,1.0],linestyle='--',color='b')
	plt.text(0.037,0.7,'Mean\n'+r'$XCO_2$'+'\nPrecision',color='b',horizontalalignment='center')
	plt.plot([0,max(bp_KBr_XCO2)],[0.97,0.97],linestyle='--',color='k')
	plt.text(0.35,0.93,'97%',color='k')
	plt.plot(bp_CaF2_XCO2,cum_func_CaF2_XCO2,color='b',label='$CaF_2$',linestyle='-',marker='None')
	plt.plot(bp_KBr_XCO2,cum_func_KBr_XCO2,color='r',label='KBr',linestyle='-',marker='None')
	plt.ylabel('Cumulative Function')
	plt.xlabel('Precision')
	plt.ylim([0.0,1.05])
	plt.legend(loc=4)
	#plt.savefig('CUM_XCO2.eps')

	plt.figure()
	plt.subplot(211)
	plt.title('$CO_2$')
	plt.plot([sigma_CO2_KBrII,sigma_CO2_KBrII],[0.0,1.0],linestyle='--',color='r')
	plt.plot([sigma_CO2_CaF2II,sigma_CO2_CaF2II],[0.0,1.0],linestyle='--',color='b')
	plt.plot(bp_CaF2_CO2,cum_func_CaF2_CO2,color='b',label='$CaF_2$',linestyle='-')
	plt.plot(bp_KBr_CO2,cum_func_KBr_CO2,color='r',label='KBr',linestyle='-')
	plt.ylabel('Cumulative Function')
	plt.ylim([0.0,1.05])
	plt.legend(loc=4)
	plt.subplot(212)
	plt.title('$O_2$')
	plt.plot([sigma_O2_KBrII,sigma_O2_KBrII],[0.0,1.0],linestyle='--',color='r')
	plt.plot([sigma_O2_CaF2II,sigma_O2_CaF2II],[0.0,1.0],linestyle='--',color='b')
	plt.plot(bp_CaF2_O2,cum_func_CaF2_O2,color='b',label='$CaF_2$',linestyle='-')
	plt.plot(bp_KBr_O2,cum_func_KBr_O2,color='r',label='KBr',linestyle='-')
	plt.ylabel('Cumulative Function')
	plt.xlabel('Standard Deviation')
	plt.ylim([0.0,1.05])
	#plt.savefig('CUM_COLS.eps')

	f_size=10

	print 'SZA MEAN',np.mean(KBr_ens['SZA']),np.mean(CaF2_ens['SZA'])
	plt.figure()
	plt.suptitle('CONTIGUOUS DETREND ('+archivo_txt+')')
	plt.subplot(311)
	plt.plot(dt_CaF2,dtr_cont_CaF2['dtr_CO2'],'.b',label='CaF2')
	plt.plot(dt_KBr,dtr_cont_KBr['dtr_CO2'],'.r',label='KBr')
	plt.ylabel('CO2 (molecules/cm2)')
	plt.xticks(rotation=45,fontsize=8)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
	plt.legend()
	plt.subplot(312)
	plt.plot(dt_CaF2,dtr_cont_CaF2['dtr_O2'],'.b',label='CaF2')
	plt.plot(dt_KBr,dtr_cont_KBr['dtr_O2'],'.r',label='KBr')
	plt.ylabel('O2 (molecules/cm2)')
	plt.xticks(rotation=45,fontsize=8)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
	plt.legend()
	plt.subplot(313)
	plt.plot(dt_CaF2,dtr_cont_CaF2['dtr_XCO2'],'.b',label='CaF2')
	plt.plot(dt_KBr,dtr_cont_KBr['dtr_XCO2'],'.r',label='KBr')
	#plt.xlabel('Epoch Time')
	plt.ylabel('XCO2 (ppm)')
	plt.xticks(rotation=45,fontsize=8)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
	#plt.legend()

	plt.figure()
	#plt.suptitle('PROFFIT ERRORS ('+archivo_txt+')')
	plt.subplot(221)
	plt.plot(dt_CaF2,CaF2_ens['LOSerrCO2'],'.b',label='CO2 CaF2')
	plt.plot(dt_KBr,KBr_ens['LOSerrCO2'],'.r',label='CO2 KBr')
	plt.ylabel('LOS ERROR')
	plt.xticks(rotation=45,fontsize=8)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
	plt.legend(fontsize=8,loc=9,ncol=2)
	plt.subplot(223)
	plt.plot(dt_CaF2,CaF2_ens['NOISEerrCO2'],'.b',label='CO2 CaF2')
	plt.plot(dt_KBr,KBr_ens['NOISEerrCO2'],'.r',label='CO2 KBr')
	plt.ylabel('NOISE ERROR')
	plt.xticks(rotation=45,fontsize=8)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
	#plt.legend(fontsize=8,loc=9,ncol=2)
	plt.subplot(222)
	plt.plot(dt_CaF2,CaF2_ens['LOSerrO2'],'.b',label='O2 CaF2')
	plt.plot(dt_KBr,KBr_ens['LOSerrO2'],'.r',label='O2 KBr')
	#plt.xlabel('Epoch Time')
	#plt.ylabel('LOS ERROR')
	plt.xticks(rotation=45,fontsize=8)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
	plt.legend(fontsize=8,loc=9,ncol=2)
	plt.subplot(224)
	plt.plot(dt_CaF2,CaF2_ens['NOISEerrO2'],'.b',label='O2 CaF2')
	plt.plot(dt_KBr,KBr_ens['NOISEerrO2'],'.r',label='O2 KBr')
	#plt.xlabel('Epoch Time')
	#plt.ylabel('NOISE ERROR')
	plt.xticks(rotation=45,fontsize=8)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %y'))
	#plt.legend(fontsize=8,loc=9,ncol=2)
	#plt.savefig('PRFERRORS_HR.eps')


##############################

COLUMNS("1211-1512_XCO2.txt")

plt.show()


