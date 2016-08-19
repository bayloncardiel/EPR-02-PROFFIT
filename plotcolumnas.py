import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.dates as mdates
import sys
import time
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'
#sys.path.append("/home/jorge/PICARRO")
#import dailymean

######################################################################
####################	TRUE SOLAR TIME
def solar_time(dattim,timezone=-6,lon=-98.6552):
	
	tst = []
	for ele in dattim:
		D = ele.timetuple().tm_yday
		gamma = ((2 * math.pi) / 365.0) * (D - 1.0 + ((ele.hour - 12.0) / 24.0))
		eqtime1 = 229.18 * (0.000075 + 0.001868 * math.cos(gamma) - 0.032077 * math.sin(gamma) \
		- 0.014615 * math.cos(2 * gamma) - 0.040849 * math.sin(2 * gamma))
		geo_time_diff =  + 4.0 * lon - 60.0 * timezone
		time_offset1 = eqtime1 + geo_time_diff
		tst1 = ele.hour * 60.0 + ele.minute + ele.second / 60.0 + time_offset1
		solar_time1 = dt.datetime.combine(ele.date(), dt.time(0)) + dt.timedelta(minutes=tst1)
		tst.append(solar_time1)
	tst = np.array(tst)
	return tst

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
####################	FILENAMES AND PATH
f_path = "/home/jorge/PROFFIT"
n_CO2 = "1211-1512_CO2err.txt"
n_O2 = "1211-1512_O2err.txt"

######################################################################
####################	READING WITH DTYPE FORMAT
dt_txt = np.dtype([('fname',np.str_,30),('epochTime',int),('SZA',float), \
('totCol',float),('RMS',float),('relRMS',float),('wnShift',float),('LOSerr',float),('NOISEerr',float)])
f_CO2 = np.genfromtxt(f_path+'/'+n_CO2,dtype = dt_txt,delimiter = " ")
f_CO2 = f_CO2[f_CO2['epochTime'].argsort()]	### SORTING BASED ON EPOCH TIME
f_O2 = np.genfromtxt(f_path+'/'+n_O2,dtype = dt_txt,delimiter = " ")
f_O2 = f_O2[f_O2['epochTime'].argsort()]	### SORTING BASED ON EPOCH TIME

dt_CO2 = np.array([dt.datetime.fromtimestamp(ele) for ele in f_CO2['epochTime']])		### DATETIME FOR CO2
dt_O2 = np.array([dt.datetime.fromtimestamp(ele) for ele in f_O2['epochTime']])			### DATETIME FOR O2

######################################################################
####################	HISTOGRAMS FOR RELATIVE RMS AND WAVENUMBER SHIFT
'''
# CO2
hist_rr_CO2,be_rr_CO2 = np.histogram(f_CO2['relRMS'],bins=50)
bp_rr_CO2 = np.array([(ele + be_rr_CO2[ii+1])/2.0 for ii,ele in enumerate(be_rr_CO2[:-1])])
hist_ws_CO2,be_ws_CO2 = np.histogram(f_CO2['wnShift'],bins=50)
bp_ws_CO2 = np.array([(ele + be_ws_CO2[ii+1])/2.0 for ii,ele in enumerate(be_ws_CO2[:-1])])

plt.figure()
plt.subplot(311)
plt.plot(bp_rr_CO2,hist_rr_CO2,'-r',label='Relative RMS CO2')
plt.ylabel('Histogram\nRelative RMS CO2')
plt.subplot(312)
plt.plot(bp_ws_CO2,hist_ws_CO2,'-b',label='Wavenumber Shift CO2')
plt.ylabel('Histogram\nWavenumber Shift CO2')
plt.subplot(313)
plt.plot(dt_CO2,f_CO2['wnShift'],'.g',label='Wavenumber Shift CO2')
plt.ylabel('Wavenumber Shift CO2')

plt.figure()
plt.suptitle('ERRORS CO2')
plt.subplot(211)
plt.plot(f_CO2['epochTime'],f_CO2['LOSerr'],'.',label='LOS')
plt.legend()
plt.subplot(212)
plt.plot(f_CO2['epochTime'],f_CO2['LOSerr'],'.',label='NOISE')
plt.legend()

# O2
hist_rr_O2,be_rr_O2 = np.histogram(f_O2['relRMS'],bins=50)
bp_rr_O2 = np.array([(ele + be_rr_O2[ii+1])/2.0 for ii,ele in enumerate(be_rr_O2[:-1])])
hist_ws_O2,be_ws_O2 = np.histogram(f_O2['wnShift'],bins=50)
bp_ws_O2 = np.array([(ele + be_ws_O2[ii+1])/2.0 for ii,ele in enumerate(be_ws_O2[:-1])])

plt.figure()
plt.subplot(311)
plt.plot(bp_rr_O2,hist_rr_O2,'-r',label='Relative RMS O2')
plt.ylabel('Histogram\nRelative RMS O2')
plt.subplot(312)
plt.plot(bp_ws_O2,hist_ws_O2,'-b',label='Wavenumber Shift O2')
plt.ylabel('Histogram\nWavenumber Shift O2')
plt.subplot(313)
plt.plot(dt_O2,f_O2['wnShift'],'.g',label='Wavenumber Shift O2')
plt.ylabel('Wavenumber Shift O2')

plt.figure()
plt.suptitle('ERRORS O2')
plt.subplot(211)
plt.plot(f_O2['epochTime'],f_O2['LOSerr'],'.',label='LOS')
plt.legend()
plt.subplot(212)
plt.plot(f_O2['epochTime'],f_O2['LOSerr'],'.',label='NOISE')
plt.legend()
'''
######################################################################
####################	FILTER CONDITIONS
'''
# 1ST FILTER CO2
print 'SHAPE CO2', f_CO2.shape
con_CO2_I = ((f_CO2['relRMS'] > 0.0) & (f_CO2['relRMS'] < 0.005)) & \
((f_CO2['wnShift'] > 3.0e-7) & (f_CO2['wnShift'] < 1.3e-6))	### ((f_CO2['relRMS'] > 0.0) & (f_CO2['relRMS'] < 0.005)) (ORIGINAL, SE CAMBIO PARA 									### 151111)
								### PRIMER FILTRO USADO ANTES DE 13/05/16
#con_CO2_I = ((f_CO2['relRMS'] > 0.0) & (f_CO2['relRMS'] < 1.0)) & \
#((f_CO2['wnShift'] > -3.0e-1) & (f_CO2['wnShift'] < 1.3e-1))	### ((f_CO2['relRMS'] > 0.0) & (f_CO2['relRMS'] < 0.005)) (ORIGINAL, SE CAMBIO PARA 
f_CO2 = f_CO2[con_CO2_I]
dt_CO2 = dt_CO2[con_CO2_I]
print 'FILTER #1 CO2', f_CO2.shape
'''
# 2ND FILTER CO2
print 'SHAPE CO2', f_CO2.shape
fac_fil_CO2 = 2.0
up_fil_CO2 = np.mean(f_CO2['totCol']) + fac_fil_CO2*np.std(f_CO2['totCol'])
dw_fil_CO2 = np.mean(f_CO2['totCol']) - fac_fil_CO2*np.std(f_CO2['totCol'])
#print 'CO2 MEAN STD',np.mean(f_CO2['totCol']),np.std(f_CO2['totCol'])
#con_CO2_II = ((f_CO2['totCol'] > dw_fil_CO2) & (f_CO2['totCol'] < up_fil_CO2))	### SEGUNDO FILTRO USADO ANTES DE 13/05/16
#f_CO2 = f_CO2[con_CO2_II]
#dt_CO2 = dt_CO2[con_CO2_II]
tsdt_CO2 = solar_time(dt_CO2)
print 'FILTER #2 CO2', f_CO2.shape
'''
# 1ST FILTER O2
print 'SHAPE O2', f_O2.shape
con_O2_I = ((f_O2['relRMS'] > 0.0) & (f_O2['relRMS'] < 0.01)) & \
((f_O2['wnShift'] > 3.0e-7) & (f_O2['wnShift'] < 1.3e-6))	### PRIMER FILTRO USADO ANTES DE 13/05/16
#con_O2_I = ((f_O2['relRMS'] > 0.0) & (f_O2['relRMS'] < 1.0)) & \
#((f_O2['wnShift'] > -3.0e-1) & (f_O2['wnShift'] < 1.3e-1))

f_O2 = f_O2[con_O2_I]
dt_O2 = dt_O2[con_O2_I]
print 'FILTER #1 O2', f_O2.shape
'''
# 2ND FILTER O2
print 'SHAPE O2', f_O2.shape
fac_fil_O2 = 2.0
up_fil_O2 = np.mean(f_O2['totCol']) + fac_fil_O2*np.std(f_O2['totCol'])
dw_fil_O2 = np.mean(f_O2['totCol']) - fac_fil_O2*np.std(f_O2['totCol'])
#print 'O2 MEAN STD',np.mean(f_O2['totCol']),np.std(f_O2['totCol'])
#con_O2_II = ((f_O2['totCol'] > dw_fil_O2) & (f_O2['totCol'] < up_fil_O2))		### SEGUNDO FILTRO USADO ANTES DE 13/05/16
#f_O2 = f_O2[con_O2_II]
#dt_O2 = dt_O2[con_O2_II]
tsdt_O2 = solar_time(dt_O2)
print 'FILTER #2 O2', f_O2.shape

######################################################################
####################	XCO2 CALCULATION
XCO2_et = np.intersect1d(f_CO2["epochTime"],f_O2["epochTime"])
XCO2_dt = np.array([dt.datetime.fromtimestamp(ele) for ele in XCO2_et])
XCO2_d = np.array([ele.date() for ele in XCO2_dt])
XCO2_tsdt = solar_time(XCO2_dt)
print 'COMMON ELEMENTS: %i' % len(XCO2_et)
print 'NUMBER OF DAYS: %i' % len(np.unique(XCO2_d))
CO2_ind = np.nonzero(np.in1d(f_CO2["epochTime"],XCO2_et))[0]
O2_ind = np.nonzero(np.in1d(f_O2["epochTime"],XCO2_et))[0]
CO2_cmm = f_CO2[CO2_ind]
O2_cmm = f_O2[O2_ind]
print 'SAME EPOCHTIME IN CO2 AND O2?',np.unique(CO2_cmm['epochTime'] == O2_cmm['epochTime'])
XCO2 = np.array([209500.0*(ele/O2_cmm['totCol'][ii]) for ii,ele in enumerate(CO2_cmm['totCol'])])

######################################################################
####################	WRITE XCO2 VALUES IN TEXT FILE
def writeXCO2txt(CO2_mat,O2_mat,XCO2):
	string_txt = ''
	for ii,ele in enumerate(XCO2):
		if CO2_mat['epochTime'][ii] == O2_mat['epochTime'][ii]:
			string_txt = string_txt+CO2_mat['fname'][ii]+' '+str(CO2_mat['epochTime'][ii])+' %.4f %.4f %.4e %.4e %.6f %.6f %.6f %.6f\n' % (CO2_mat['SZA'][ii],ele,CO2_mat['totCol'][ii],O2_mat['totCol'][ii],CO2_mat['LOSerr'][ii],CO2_mat['NOISEerr'][ii],O2_mat['LOSerr'][ii],O2_mat['NOISEerr'][ii])
	min_date = dt.datetime.fromtimestamp(min(CO2_mat['epochTime']))
	max_date = dt.datetime.fromtimestamp(max(CO2_mat['epochTime']))
	file_XCO2 = open('%02i%02i-%02i%02i_XCO2.txt' % (min_date.year-2000,min_date.month,max_date.year-2000,max_date.month),'w')
	file_XCO2.write(string_txt)
	file_XCO2.close()
	print 'ARCHIVO %02i%02i-%02i%02i_XCO2.txt GENERADO' % (min_date.year-2000,min_date.month,max_date.year-2000,max_date.month)
writeXCO2txt(CO2_cmm,O2_cmm,XCO2)

######################################################################
####################	DAILY AND MONTLY MEANS
# D M CO2
d_CO2 = np.array([ele.date() for ele in dt_CO2])
unq_d_CO2 = np.unique(d_CO2)
dailymean_CO2 = np.array([np.mean(f_CO2['totCol'][d_CO2 == ele]) for ele in unq_d_CO2])
dailystd_CO2 = np.array([np.std(f_CO2['totCol'][d_CO2 == ele])/np.sqrt(len(f_CO2['totCol'][d_CO2 == ele])) for ele in unq_d_CO2])
# M M CO2
my_CO2 = np.array(['%04i-%02i' % (ele.year,ele.month) for ele in d_CO2])
unq_my_CO2 = np.unique(my_CO2)
monthlymean_CO2 = np.array([np.mean(f_CO2['totCol'][my_CO2 == ele]) for ele in unq_my_CO2])
monthlystd_CO2 = np.array([np.std(f_CO2['totCol'][my_CO2 == ele])/np.sqrt(len(f_CO2['totCol'][my_CO2 == ele])) for ele in unq_my_CO2])
#dtrend_CO2 = np.array([f_CO2['totCol'][ii] - dailymean_CO2[unq_d_CO2 == ele] for ii,ele in enumerate(d_CO2)])		### DAILY MEAN DETRENDING
dtrend_CO2 = np.array([f_CO2['totCol'][ii] - monthlymean_CO2[unq_my_CO2 == ele] for ii,ele in enumerate(my_CO2)])	### MONTHLY MEAN DETRENDING
unq_my_CO2 = np.array([dt.datetime.strptime(ele+'-15','%Y-%m-%d') for ele in unq_my_CO2])

# D M O2
d_O2 = np.array([ele.date() for ele in dt_O2])
unq_d_O2 = np.unique(d_O2)
dailymean_O2 = np.array([np.mean(f_O2['totCol'][d_O2 == ele]) for ele in unq_d_O2])
dailystd_O2 = np.array([np.std(f_O2['totCol'][d_O2 == ele])/np.sqrt(len(f_O2['totCol'][d_O2 == ele])) for ele in unq_d_O2])
# M M O2
my_O2 = np.array(['%04i-%02i' % (ele.year,ele.month) for ele in d_O2])
unq_my_O2 = np.unique(my_O2)
monthlymean_O2 = np.array([np.mean(f_O2['totCol'][my_O2 == ele]) for ele in unq_my_O2])
monthlystd_O2 = np.array([np.std(f_O2['totCol'][my_O2 == ele])/np.sqrt(len(f_O2['totCol'][my_O2 == ele])) for ele in unq_my_O2])
#dtrend_O2 = np.array([f_O2['totCol'][ii] - dailymean_O2[unq_d_O2 == ele] for ii,ele in enumerate(d_O2)])	### DAILY MEAN DETRENDING
dtrend_O2 = np.array([f_O2['totCol'][ii] - monthlymean_O2[unq_my_O2 == ele] for ii,ele in enumerate(my_O2)])	### MONTHLY MEAN DETRENDING
unq_my_O2 = np.array([dt.datetime.strptime(ele+'-15','%Y-%m-%d') for ele in unq_my_O2])

# D M XCO2
XCO2_d = np.array([ele.date() for ele in XCO2_dt])
XCO2_unq_d = np.unique(XCO2_d)
XCO2_dailymean = np.array([np.mean(XCO2[XCO2_d == ele]) for ele in XCO2_unq_d])
XCO2_dailystd = np.array([np.std(XCO2[XCO2_d == ele])/np.sqrt(len(XCO2[XCO2_d == ele])) for ele in XCO2_unq_d])
# M M XCO2
XCO2_my = np.array(['%04i-%02i' % (ele.year,ele.month) for ele in XCO2_dt])
XCO2_unq_my = np.unique(XCO2_my)
XCO2_monthlymean = np.array([np.mean(XCO2[XCO2_my == ele]) for ele in XCO2_unq_my])
XCO2_monthlystd = np.array([np.std(XCO2[XCO2_my == ele])/np.sqrt(len(XCO2[XCO2_my == ele])) for ele in XCO2_unq_my])
XCO2_dtrend = np.array([XCO2[ii] - XCO2_dailymean[XCO2_unq_d == ele] for ii,ele in enumerate(XCO2_d)])		### DAILY MEAN DETRENDING
#XCO2_dtrend = np.array([XCO2[ii] - XCO2_monthlymean[XCO2_unq_my == ele] for ii,ele in enumerate(XCO2_my)])	### MONTHLY MEAN DETRENDING
XCO2_unq_my = np.array([dt.datetime.strptime(ele+'-15','%Y-%m-%d') for ele in XCO2_unq_my])

######################################################################
####################	HOURLY AND WEEKDAY MEANS
# H M CO2
hr_CO2 = np.array([ele.hour for ele in dt_CO2])
tshr_CO2 = np.array([ele.hour for ele in tsdt_CO2])
unq_hr_CO2 = np.unique(hr_CO2)
unq_tshr_CO2 = np.unique(tshr_CO2)
hourlymean_CO2 = np.array([np.mean(f_CO2['totCol'][hr_CO2 == ele]) for ele in unq_hr_CO2])
hourlystd_CO2 = np.array([np.std(f_CO2['totCol'][hr_CO2 == ele])/np.sqrt(len(f_CO2['totCol'][hr_CO2 == ele])) for ele in unq_hr_CO2])
tshourlymean_CO2 = np.array([np.mean(f_CO2['totCol'][tshr_CO2 == ele]) for ele in unq_tshr_CO2])
tshourlystd_CO2 = np.array([np.std(f_CO2['totCol'][tshr_CO2 == ele])/np.sqrt(len(f_CO2['totCol'][tshr_CO2 == ele])) for ele in unq_tshr_CO2])
# WD M CO2
wd_CO2 = np.array([ele.weekday() for ele in dt_CO2])
unq_wd_CO2 = np.unique(wd_CO2)
weekdaymean_CO2 = np.array([np.mean(f_CO2['totCol'][wd_CO2 == ele]) for ele in unq_wd_CO2])
weekdaystd_CO2 = np.array([np.std(f_CO2['totCol'][wd_CO2 == ele])/np.sqrt(len(f_CO2['totCol'][wd_CO2 == ele])) for ele in unq_wd_CO2])

# H M O2
hr_O2 = np.array([ele.hour for ele in dt_O2])
tshr_O2 = np.array([ele.hour for ele in tsdt_O2])
unq_hr_O2 = np.unique(hr_O2)
unq_tshr_O2 = np.unique(tshr_O2)
hourlymean_O2 = np.array([np.mean(f_O2['totCol'][hr_O2 == ele]) for ele in unq_hr_O2])
hourlystd_O2 = np.array([np.std(f_O2['totCol'][hr_O2 == ele])/np.sqrt(len(f_O2['totCol'][hr_O2 == ele])) for ele in unq_hr_O2])
tshourlymean_O2 = np.array([np.mean(f_O2['totCol'][tshr_O2 == ele]) for ele in unq_tshr_O2])
tshourlystd_O2 = np.array([np.std(f_O2['totCol'][tshr_O2 == ele])/np.sqrt(len(f_O2['totCol'][tshr_O2 == ele])) for ele in unq_tshr_O2])
# WD M O2
wd_O2 = np.array([ele.weekday() for ele in dt_O2])
unq_wd_O2 = np.unique(wd_O2)
weekdaymean_O2 = np.array([np.mean(f_O2['totCol'][wd_O2 == ele]) for ele in unq_wd_O2])
weekdaystd_O2 = np.array([np.std(f_O2['totCol'][wd_O2 == ele])/np.sqrt(len(f_O2['totCol'][wd_O2 == ele])) for ele in unq_wd_O2])

# H M XCO2
XCO2_hr = np.array([ele.hour for ele in XCO2_dt])
XCO2_tshr = np.array([ele.hour for ele in XCO2_tsdt])
XCO2_unq_hr = np.unique(XCO2_hr)
XCO2_unq_tshr = np.unique(XCO2_tshr)
XCO2_hourlymean = np.array([np.mean(XCO2[XCO2_hr == ele]) for ele in XCO2_unq_hr])
XCO2_hourlystd = np.array([np.std(XCO2[XCO2_hr == ele])/np.sqrt(len(XCO2[XCO2_hr == ele])) for ele in XCO2_unq_hr])
XCO2_tshourlymean = np.array([np.mean(XCO2[XCO2_tshr == ele]) for ele in XCO2_unq_tshr])
XCO2_tshourlystd = np.array([np.std(XCO2[XCO2_tshr == ele])/np.sqrt(len(XCO2[XCO2_tshr == ele])) for ele in XCO2_unq_tshr])
# WD M XCO2
XCO2_wd = np.array([ele.weekday() for ele in XCO2_dt])
XCO2_unq_wd = np.unique(XCO2_wd)
XCO2_weekdaymean = np.array([np.mean(XCO2[XCO2_wd == ele]) for ele in XCO2_unq_wd])
XCO2_weekdaystd = np.array([np.std(XCO2[XCO2_wd == ele])/np.sqrt(len(XCO2[XCO2_wd == ele])) for ele in XCO2_unq_wd])

######################################################################
####################	FOURIER SERIES
# XCO2
#min_d = np.amin(XCO2_unq_d) + dt.timedelta(days = (max(XCO2_unq_d) - min(XCO2_unq_d)).days/2)
#diff_d = np.array([(ele - min_d).days for ele in XCO2_unq_d])
#diff_d = np.array(decimal_year(XCO2_unq_d))
#numday_d = np.array([ele.strftime('%j') for ele in XCO2_unq_d])
#yr_XCO2_dailymean = XCO2_dailymean[(diff_d >= 2015.0) & (diff_d < 2016.0)]
#yr_numday_d = numday_d[(diff_d >= 2015.0) & (diff_d < 2016.0)]
#print yr_numday_d[yr_XCO2_dailymean == min(yr_XCO2_dailymean)]
#print yr_numday_d[yr_XCO2_dailymean == max(yr_XCO2_dailymean)]
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
	print 'For %s (%i): a0 = %.6e, alpha = %.6e, a1 = %.6e, a2 = %.6e, b1 = %.6e, b2 = %.6e' % (species,len(meas_val),coef[0],coef[1],coef[2],coef[4],coef[3],coef[5])
	#fourier_fun = np.zeros(len(XCO2_unq_d))
	fourier_fun = np.zeros(len(diff_d))
	trend = np.array([coef[0] + coef[1]*ele for ele in diff_d])
	for ii,ele in enumerate(K[:,0]):
		for jj,cof in enumerate(coef):
			fourier_fun[ii] = fourier_fun[ii]+K[ii,jj]*cof
	return fourier_fun,trend,coef
datetime_CO2 = [dt.datetime.fromtimestamp(ala) for ala in CO2_cmm['epochTime']]
datetime_O2 = [dt.datetime.fromtimestamp(ele) for ele in O2_cmm['epochTime']]
diff_d_CO2 = np.array(decimal_year(datetime_CO2))
diff_d_O2 = np.array(decimal_year(datetime_O2))
diff_d_XCO2 = np.array(decimal_year(XCO2_dt))
fourier_fun_CO2,trend_CO2,coef_CO2 = fourier_series(diff_d_CO2,CO2_cmm['totCol'],'CO2')
fourier_fun_O2,trend_O2,coef_O2 = fourier_series(diff_d_O2,O2_cmm['totCol'],'O2')
fourier_fun_XCO2,trend_XCO2,coef_XCO2 = fourier_series(diff_d_XCO2,XCO2,'XCO2')
dtr_ff_CO2 = np.array([CO2_cmm['totCol'][aa]-ala for aa,ala in enumerate(fourier_fun_CO2)])
dtr_ff_O2 = np.array([O2_cmm['totCol'][ee]-ele for ee,ele in enumerate(fourier_fun_O2)])
dtr_ff_XCO2 = np.array([XCO2[ii]-ili for ii,ili in enumerate(fourier_fun_XCO2)])
abspath_fname = []
for ee,ele in enumerate(CO2_cmm['fname']):
	if dt.datetime.fromtimestamp(CO2_cmm['epochTime'][ee]).date() <= dt.date(year=2015,month=9,day=24):
		abspath_fname.append('/home/D2_PROFFIT/Jorge_Files/O2_v1/%s/%s' % (ele.split('.')[0],ele))
	else:
		abspath_fname.append('/home/D2_RESULTS/PROFFIT_results/ALTZ/O2_v1/%s/%s' % (ele.split('.')[0],ele))
abspath_fname = np.array(abspath_fname)
dset_dtr_XCO2 = np.empty(len(CO2_cmm['epochTime']),dtype=[('filenamehdf',np.str_,100),('tepoch',int),('totcol',float)])
dset_dtr_XCO2['filenamehdf'] = abspath_fname
dset_dtr_XCO2['tepoch'] = CO2_cmm['epochTime']
o2_mean = np.mean(O2_cmm['totCol'])
#dset_dtr_XCO2['totcol'] = (dtr_ff_XCO2)*(o2_mean/209500.0)
dset_dtr_XCO2['totcol'] = (dtr_ff_O2)
#np.save('/home/D2_RESULTS/PROFFIT_results/TIMESERIES_analysis/npdata/dtrO2_v1_QA.npy',dset_dtr_XCO2)
#mat_dtr = np.column_stack((CO2_cmm['fname'],CO2_cmm['epochTime'],CO2_cmm['SZA'],XCO2,CO2_cmm['totCol'],O2_cmm['totCol'],dtr_ff_XCO2,dtr_ff_CO2,dtr_ff_O2))
#np.savetxt('dtrcol.txt',mat_dtr,fmt="%s")		### MATRIZ CON XCO2, CO2 Y O2 DETRENDED


######################################################################
####################	PLOTS

f_size = 12
plt_color = "#FF4F2F"

plt.figure()
plt.plot(XCO2_dt,XCO2,color=plt_color,marker='.',linestyle='None',label=None)
plt.plot(XCO2_dt,fourier_fun_XCO2,linestyle='-',marker='None',linewidth=3,color='k',label=None)
plt.plot(XCO2_dt,trend_XCO2,linestyle='--',marker='None',linewidth=1,color='k',label='%.2f ppm/yr' % (coef_XCO2[1]))
#plt.errorbar(XCO2_unq_d,XCO2_dailymean,yerr=XCO2_dailystd,color=plt_color,marker='.',linestyle='None',label=None)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))	#[1,15]
#plt.ylim([390.0,405.0])
plt.yticks(fontsize=f_size)
plt.xticks(rotation=45,fontsize=f_size)		
plt.ylabel('$XCO_2\/(ppm)$',fontsize=f_size)
#plt.legend(ncol=3,fontsize=8,loc=9)
#plt.savefig('ANNUALCYCLE_XCO2III.jpg')

plt.figure()
ax1 = plt.subplot(211)
plt.plot(datetime_CO2,CO2_cmm['totCol'],color=plt_color,marker='.',linestyle='None',label=None)
plt.plot(datetime_CO2,fourier_fun_CO2,linestyle='-',marker='None',linewidth=3,color='k',label=None)
plt.plot(datetime_CO2,trend_CO2,linestyle='--',marker='None',linewidth=1,color='k',label='%.2f ppm/yr' % (coef_CO2[1]))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))	#[1,15]
plt.yticks(fontsize=f_size)
plt.xticks(rotation=45,fontsize=f_size)		
plt.ylabel('$CO_2\/(mol/cm^2)$',fontsize=f_size)
ax2 = plt.subplot(212,sharex=ax1)
plt.plot(datetime_O2,O2_cmm['totCol'],color=plt_color,marker='.',linestyle='None',label=None)
plt.plot(datetime_O2,fourier_fun_O2,linestyle='-',marker='None',linewidth=3,color='k',label=None)
plt.plot(datetime_O2,trend_O2,linestyle='--',marker='None',linewidth=1,color='k',label='%.2f ppm/yr' % (coef_O2[1]))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))	#[1,15]
plt.yticks(fontsize=f_size)
plt.xticks(rotation=45,fontsize=f_size)		
plt.ylabel('$O_2\/(mol/cm^2)$',fontsize=f_size)

plt.figure()
ax3 = plt.subplot(311)
plt.plot(datetime_CO2,dtr_ff_CO2,color=plt_color,marker='.',linestyle='None',label=None)
plt.ylabel('$CO_2\/Detrended$')
ax4 = plt.subplot(312,sharex=ax3)
plt.plot(datetime_O2,dtr_ff_O2,color=plt_color,marker='.',linestyle='None',label=None)
plt.ylabel('$O_2\/Detrended$')
ax5 = plt.subplot(313,sharex=ax3)
plt.plot(XCO2_dt,dtr_ff_XCO2,color=plt_color,marker='.',linestyle='None',label=None)
plt.ylabel('$XCO_2\/Detrended$')

print 'MIN (%.4f) MAX (%.4f) DIFF WEEKLY CYCLE: %.4f' % (min(XCO2_weekdaymean),max(XCO2_weekdaymean),max(XCO2_weekdaymean)-min(XCO2_weekdaymean))
print 'MIN (%.4f) MAX (%.4f) DIFF DIURNAL CYCLE: %.4f' % (min(XCO2_tshourlymean[1:-1]),max(XCO2_tshourlymean[1:-1]),max(XCO2_tshourlymean[1:-1])-min(XCO2_tshourlymean[1:-1]))
'''
plt.figure()
#plt.plot(XCO2_unq_d,fourier_fun,'--k')
plt.errorbar(unq_d_CO2,dailymean_CO2/(10**21),yerr=dailystd_CO2/(10**21),color=plt_color,marker='.',linestyle='None',label='DAILY MEAN CO2')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))	#[1,15]
#plt.ylim([390.0,405.0])
plt.yticks(fontsize=f_size)
plt.xticks(rotation=45,fontsize=f_size)		
plt.ylabel('$CO_2\/(10^{21}\/molecules\/cm^{-1}$)',fontsize=f_size)
#plt.savefig('ANNUALCYCLE_CO2.jpg')

plt.figure()
#plt.plot(XCO2_unq_d,fourier_fun,'--k')
plt.errorbar(unq_d_O2,dailymean_O2/(10**24),yerr=dailystd_O2/(10**24),color=plt_color,marker='.',linestyle='None',label='DAILY MEAN O2')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))	#[1,15]
#plt.ylim([390.0,405.0])
plt.yticks(fontsize=f_size)
plt.xticks(rotation=45,fontsize=f_size)		
plt.ylabel('$O_2\/(10^{24}\/molecules\/cm^{-1}$)',fontsize=f_size)
#plt.savefig('ANNUALCYCLE_O2.jpg')
'''

plt.figure()
plt.errorbar(XCO2_unq_wd,XCO2_weekdaymean,yerr=XCO2_weekdaystd,color=plt_color,marker='.',linestyle='-',label='WEEKDAY MEAN XCO2')
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
plt.yticks(fontsize=f_size)
plt.xticks(XCO2_unq_wd,['Mon','Tue','Wed','Thu','Fri','Sat','Sun'],fontsize=f_size)#,rotation=45)
plt.ylabel('$XCO_2\/(ppm)$',fontsize=f_size)
plt.xlim([min(XCO2_unq_wd)-0.2,max(XCO2_unq_wd)+0.2])
plt.ylim([394.0,397.0])
plt.savefig('WEEKCYCLE_XCO2.eps')
'''
plt.figure()
plt.errorbar(unq_wd_CO2,weekdaymean_CO2/(10**21),yerr=weekdaystd_CO2/(10**21),color=plt_color,marker='.',linestyle='-',label='WEEKDAY MEAN CO2')
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
plt.yticks(fontsize=f_size)
plt.xticks(unq_wd_CO2,['Mon','Tue','Wed','Thu','Fri','Sat','Sun'],rotation=45,fontsize=f_size)
plt.ylabel('$CO_2\/(10^{21}\/molecules\/cm^{-1})$',fontsize=f_size)
plt.xlim([min(unq_wd_CO2)-0.2,max(unq_wd_CO2)+0.2])
plt.ylim([5.450,5.480])
plt.savefig('WEEKCYCLE_CO2.jpg')

plt.figure()
plt.errorbar(unq_wd_O2,weekdaymean_O2/(10**24),yerr=weekdaystd_O2/(10**24),color=plt_color,marker='.',linestyle='-',label='WEEKDAY MEAN O2')
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
plt.yticks(fontsize=f_size)
plt.xticks(unq_wd_O2,['Mon','Tue','Wed','Thu','Fri','Sat','Sun'],rotation=45,fontsize=f_size)
plt.ylabel('$O_2\/(10^{24}\/molecules\/cm^{-1})$',fontsize=f_size)
plt.xlim([min(unq_wd_CO2)-0.2,max(unq_wd_CO2)+0.2])
plt.ylim([2.890,2.902])
plt.savefig('WEEKCYCLE_O2.jpg')
'''
plt.figure()
#plt.errorbar(XCO2_unq_hr,XCO2_hourlymean,yerr=XCO2_hourlystd,color=plt_color,marker='.',linestyle='-',label='HOURLY MEAN XCO2')
plt.errorbar(XCO2_unq_tshr[1:-1],XCO2_tshourlymean[1:-1],yerr=XCO2_tshourlystd[1:-1],color=plt_color,marker='.',linestyle='-',label='HOURLY MEAN XCO2 TS')
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
plt.yticks(fontsize=f_size)
plt.xticks(XCO2_unq_tshr[1:-1],["%i:00" % ele for ele in XCO2_unq_tshr[1:-1]],fontsize=f_size)#,rotation=45)
plt.xlabel('$True\/Solar\/Time$',fontsize=f_size)
plt.ylabel('$XCO_2\/(ppm)$',fontsize=f_size)
plt.xlim([min(XCO2_unq_tshr[1:-1])-0.2,max(XCO2_unq_tshr[1:-1])+0.2])
plt.ylim([min(XCO2_tshourlymean[1:-1])*0.999,max(XCO2_tshourlymean[1:-1])*1.002]) #plt.ylim([394.0,397.0])
plt.savefig('DIURNALCYCLETST_XCO2.eps')

#plt.figure()
#plt.plot(yr_numday_d,yr_XCO2_dailymean,linestyle='None',marker='*',linewidth=3,color='k',label=None)

plt.show()
'''
plt.figure()
#plt.errorbar(unq_hr_CO2,hourlymean_CO2/(10**21),yerr=hourlystd_CO2/(10**21),color=plt_color,marker='.',linestyle='-',label='HOURLY MEAN CO2')
plt.errorbar(unq_tshr_CO2,tshourlymean_CO2/(10**21),yerr=tshourlystd_CO2/(10**21),color=plt_color,marker='.',linestyle='-',label='HOURLY MEAN CO2 TS')
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
plt.yticks(fontsize=f_size)
plt.xticks(unq_tshr_CO2,["%i:00" % ele for ele in unq_tshr_CO2],rotation=45,fontsize=f_size)
plt.xlabel('$True\/Solar\/Time$',fontsize=f_size)
plt.ylabel('$CO_2\/(10^{21}\/molecules\/cm^{-1})$',fontsize=f_size)
plt.xlim([min(unq_tshr_CO2)-0.2,max(unq_tshr_CO2)+0.2])
plt.ylim([(min(tshourlymean_CO2)/(10**21))*0.999,(max(tshourlymean_CO2)/(10**21))*1.001])#plt.ylim([5.450,5.480])
#plt.savefig('DIURNALCYCLETST_CO2.eps')

plt.figure()
#plt.errorbar(unq_hr_O2,hourlymean_O2/(10**24),yerr=hourlystd_O2/(10**24),color=plt_color,marker='.',linestyle='-',label='HOURLY MEAN O2')
plt.errorbar(unq_tshr_O2,tshourlymean_O2/(10**24),yerr=tshourlystd_O2/(10**24),color=plt_color,marker='.',linestyle='-',label='HOURLY MEAN O2 TS')
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
plt.yticks(fontsize=f_size)
plt.xticks(unq_tshr_O2,["%i:00" % ele for ele in unq_tshr_O2],rotation=45,fontsize=f_size)
plt.xlabel('$True\/Solar\/Time$',fontsize=f_size)
plt.ylabel('$O_2\/(10^{24}\/molecules\/cm^{-1})$',fontsize=f_size)
plt.xlim([min(unq_tshr_O2)-0.2,max(unq_tshr_O2)+0.2])
plt.ylim([(min(tshourlymean_O2)/(10**24))*0.999,(max(tshourlymean_O2)/(10**24))*1.003])#plt.ylim([2.890,2.902])
#plt.savefig('DIURNALCYCLETST_O2.eps')
'''
'''
plt.figure()
plt.suptitle('2015-03-06')
plt.subplot(311)
plt.plot(dt_CO2[d_CO2 == dt.date(2015,3,06)],f_CO2['totCol'][d_CO2 == dt.date(2015,3,06)],'.r',label='CO2')
plt.ylabel('CO2\n(mol/cm2)')
plt.subplot(312)
plt.plot(dt_O2[d_O2 == dt.date(2015,3,06)],f_O2['totCol'][d_O2 == dt.date(2015,3,06)],'.r',label='O2')
plt.ylabel('O2\n(mol/cm2)')
plt.subplot(313)
plt.plot(XCO2_dt[XCO2_d == dt.date(2015,3,06)],XCO2[XCO2_d == dt.date(2015,3,06)],'.r',label='XCO2')
plt.ylabel('XCO2\n(ppm)')

plt.figure()
plt.suptitle('2015-03-06')
plt.subplot(311)
plt.plot(dt_CO2,f_CO2['totCol'],'.r',label='CO2')
plt.ylabel('CO2\n(mol/cm2)')
plt.subplot(312)
plt.plot(dt_O2,f_O2['totCol'],'.r',label='O2')
plt.ylabel('O2\n(mol/cm2)')
plt.subplot(313)
plt.plot(XCO2_dt,XCO2,'.r',label='XCO2')
plt.ylabel('XCO2\n(ppm)')

plt.show()
#plt.plot(dt_CO2,dtrend_CO2,'.r',label='DETRENDED CO2')
#plt.ylabel('Detrended CO2\n(mol/cm2)')
#plt.subplot(312)
#plt.errorbar(unq_d_CO2,dailymean_CO2,yerr=dailystd_CO2,color='b',marker='.',linestyle='None',label='DAILY MEAN CO2')
#plt.ylabel('Daily Mean CO2\n(mol/cm2)')
#plt.errorbar(unq_hr_CO2,hourlymean_CO2,yerr=hourlystd_CO2,color='b',marker='.',linestyle='--',label='HOURLY MEAN CO2')
#plt.ylabel('Hourly Mean CO2\n(mol/cm2)')
#plt.subplot(313)
#plt.errorbar(unq_my_CO2,monthlymean_CO2,yerr=monthlystd_CO2,color='g',marker='.',linestyle='--',label='MONTHLY MEAN CO2')
#plt.ylabel('Monthly Mean CO2\n(mol/cm2)')
#plt.errorbar(unq_wd_CO2,weekdaymean_CO2,yerr=weekdaystd_CO2,color='g',marker='.',linestyle='--',label='WEEKDAY MEAN CO2')
#plt.ylabel('Weekday Mean CO2\n(mol/cm2)')

plt.figure()
plt.suptitle('O2 TOTAL COLUMN')
#plt.subplot(311)
plt.plot(dt_O2,f_O2['totCol'],'.r',label='O2')
plt.ylabel('O2\n(mol/cm2)')
#plt.plot(dt_O2,dtrend_O2,'.r',label='DETRENDED O2')
#plt.ylabel('Detrended O2\n(mol/cm2)')
#plt.subplot(312)
#plt.errorbar(unq_d_O2,dailymean_O2,yerr=dailystd_O2,color='b',marker='.',linestyle='None',label='DAILY MEAN O2')
#plt.ylabel('Daily Mean O2\n(mol/cm2)')
#plt.errorbar(unq_hr_O2,hourlymean_O2,yerr=hourlystd_O2,color='b',marker='.',linestyle='--',label='HOURLY MEAN O2')
#plt.ylabel('Hourly Mean O2\n(mol/cm2)')
#plt.subplot(313)
#plt.errorbar(unq_my_O2,monthlymean_O2,yerr=monthlystd_O2,color='g',marker='.',linestyle='--',label='MONTHLY MEAN O2')
#plt.ylabel('Monthly Mean O2\n(mol/cm2)')
#plt.errorbar(unq_wd_O2,weekdaymean_O2,yerr=weekdaystd_O2,color='g',marker='.',linestyle='--',label='WEEKDAY MEAN O2')
#plt.ylabel('Weekday Mean O2\n(mol/cm2)')

plt.figure()
plt.suptitle('XCO2')
#plt.subplot(311)
plt.plot(XCO2_dt,XCO2,'.r',label='XCO2')
#plt.plot(XCO2_unq_d,fourier_fun,'--k')
plt.ylabel('XCO2\n(ppm)')
#plt.plot(XCO2_dt,XCO2_dtrend,'.r',label='DETRENDED XCO2')
#plt.ylabel('Detrended XCO2\n(ppm)')
#plt.subplot(312)
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
#plt.errorbar(XCO2_unq_d,XCO2_dailymean,yerr=XCO2_dailystd,color='b',marker='.',linestyle='None',label='DAILY MEAN XCO2')
#plt.xticks(rotation=45,fontsize=8)		
#plt.ylabel('Daily Mean XCO2\n(ppm)')
#plt.errorbar(XCO2_unq_hr,XCO2_hourlymean,yerr=XCO2_hourlystd,color='b',marker='.',linestyle='--',label='HOURLY MEAN XCO2')
#plt.ylabel('Hourly Mean XCO2\n(ppm)')
#plt.subplot(313)
#plt.errorbar(XCO2_unq_my,XCO2_monthlymean,yerr=XCO2_monthlystd,color='g',marker='.',linestyle='--',label='MONTHLY MEAN XCO2')
#plt.ylabel('Monthly Mean XCO2\n(ppm)')
#plt.errorbar(XCO2_unq_wd,XCO2_weekdaymean,yerr=XCO2_weekdaystd,color='g',marker='.',linestyle='--',label='WEEKDAY MEAN XCO2')
#plt.ylabel('Weekday Mean XCO2\n(ppm)')
'''
