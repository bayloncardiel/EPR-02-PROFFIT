import numpy as np
import glob
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.dates as mdates
rcParams['mathtext.default'] = 'regular'


p_path = '/home/operador/Desktop/jorge/PROFFIT/PRESION_CCA/'

list_p_files = glob.glob(p_path+'*UNAM*.dat')
list_p_files.sort()	
p_datetime = []
p_pressure = []
for ele in list_p_files:
	f = open(ele,'r')
	p_lines = f.readlines()
	f.close()
	
	for ili in p_lines[4:]:
		inf_ch = ili.split(',')
		p_dt = dt.datetime.strptime(inf_ch[0].split('"')[1],"%Y-%m-%d %H:%M:%S")
		p_press = float(inf_ch[14])
		if p_press > 510.0:
			p_datetime.append(p_dt)
			p_pressure.append(p_press)



p_datetime = np.array(p_datetime)
idx_p = np.argsort(p_datetime)
p_pressure = np.array(p_pressure)
p_datetime = p_datetime[idx_p]
p_pressure = p_pressure[idx_p]
p_date = np.array([tlt.date() for tlt in p_datetime])
p_hour = np.array([wlw.time().hour for wlw in p_datetime])

CO2_path = '/home/operador/Desktop/jorge/PROFFIT/'
CO2_file = '1603-1606_CO2_EM27.txt'
dt_txt = np.dtype([('fname',np.str_,30),('epochTime',int),('SZA',float), \
('totCol',float),('RMS',float),('relRMS',float),('wnShift',float),('LOSerr',float),('NOISEerr',float)])
CO2_mat = np.genfromtxt(CO2_path+CO2_file,dtype=dt_txt)

O2_path = '/home/operador/Desktop/jorge/PROFFIT/'
O2_file = '1603-1606_O2_EM27.txt'
dt_txt = np.dtype([('fname',np.str_,30),('epochTime',int),('SZA',float), \
('totCol',float),('RMS',float),('relRMS',float),('wnShift',float),('LOSerr',float),('NOISEerr',float)])
O2_mat = np.genfromtxt(O2_path+O2_file,dtype=dt_txt)

### FILTERS

print 'NO FILTER CO2: %i, O2: %i' % (len(CO2_mat),len(O2_mat))
'''
con_CO2_I = ((CO2_mat['relRMS'] > 0.0) & (CO2_mat['relRMS'] < 0.003))
con_O2_I = ((O2_mat['relRMS'] > 0.0) & (O2_mat['relRMS'] < 0.003))
CO2_mat = CO2_mat[con_CO2_I]
O2_mat = O2_mat[con_O2_I]

print 'AFTER FILTER I CO2: %i, O2: %i' % (len(CO2_mat),len(O2_mat))

mn_CO2 = np.mean(CO2_mat['totCol'])
std_CO2 = np.std(CO2_mat['totCol'])
con_CO2_II = ((CO2_mat['totCol'] > mn_CO2-2*std_CO2) & (CO2_mat['totCol'] < mn_CO2+2*std_CO2))
mn_O2 = np.mean(O2_mat['totCol'])
std_O2 = np.std(O2_mat['totCol'])
con_O2_II = ((O2_mat['totCol'] > mn_O2-2*std_O2) & (O2_mat['totCol'] < mn_O2+2*std_O2))
CO2_mat = CO2_mat[con_CO2_II]
O2_mat = O2_mat[con_O2_II]

print 'AFTER FILTER II CO2: %i, O2: %i' % (len(CO2_mat),len(O2_mat))
'''
###
### XCO2

XCO2_et = np.intersect1d(CO2_mat["epochTime"],O2_mat["epochTime"])
XCO2_dt = np.array([dt.datetime.fromtimestamp(ele) for ele in XCO2_et])
XCO2_d = np.array([ele.date() for ele in XCO2_dt])
print 'COMMON ELEMENTS: %i' % len(XCO2_et)
print 'NUMBER OF DAYS: %i' % len(np.unique(XCO2_d))
CO2_ind = np.nonzero(np.in1d(CO2_mat["epochTime"],XCO2_et))[0]
O2_ind = np.nonzero(np.in1d(O2_mat["epochTime"],XCO2_et))[0]
CO2_mat = CO2_mat[CO2_ind]
O2_mat = O2_mat[O2_ind]
print 'SAME EPOCHTIME IN CO2 AND O2?',np.unique(CO2_mat['epochTime'] == O2_mat['epochTime'])
XCO2 = np.array([209500.0*(ele/O2_mat['totCol'][ii]) for ii,ele in enumerate(CO2_mat['totCol'])])

###

CO2_datetime = np.array([dt.datetime.fromtimestamp(qlq) for qlq in CO2_mat['epochTime']])
CO2_date = np.array([klk.date() for klk in CO2_datetime])
CO2_col = CO2_mat['totCol']

O2_datetime = np.array([dt.datetime.fromtimestamp(qlq) for qlq in O2_mat['epochTime']])
O2_date = np.array([klk.date() for klk in O2_datetime])
O2_col = O2_mat['totCol']

print 'Dates: Min: %s, Max: %s'%(min(O2_datetime),max(O2_datetime))
print 'CO2: Min: %s, Max: %s'%(min(CO2_col),max(CO2_col))
print 'O2: Min: %s, Max: %s'%(min(O2_col),max(O2_col))
print 'XCO2: Min: %s, Max: %s'%(min(XCO2),max(XCO2))

### BOOLEAN CONDITIONS ON O2 COLUMNS AND PRESSURE
'''
fil_date = dt.date(2016,05,06)

bool_date_p = p_date == fil_date
bool_hour_p = (p_hour >= (min(O2_datetime).time().hour - 0)) & (p_hour <= (max(O2_datetime).time().hour + 1))

p_datetime = p_datetime[bool_date_p]
p_hour = p_hour[bool_date_p] 
p_pressure = p_pressure[bool_date_p]

fil_p_datetime = p_datetime[bool_hour_p]
fil_p_pressure = p_pressure[bool_hour_p]

bool_date_CO2 = CO2_date == fil_date
bool_date_O2 = O2_date == fil_date

fil_CO2_datetime = CO2_datetime[bool_date_CO2]
fil_CO2_col = CO2_col[bool_date_CO2]

fil_O2_datetime = O2_datetime[bool_date_O2]
fil_O2_col = O2_col[bool_date_O2]
'''
###
f_size=10

plt.figure(1)
ax1 = plt.subplot(311)
plt.plot(XCO2_dt,XCO2,linestyle='None',marker='.',color='r')
plt.xticks(fontsize=f_size)
plt.yticks(fontsize=f_size)		
plt.ylabel('$XCO_2$',fontsize=f_size)
ax2 = plt.subplot(312,sharex=ax1)
plt.plot(CO2_datetime,CO2_col/10**21,linestyle='None',marker='.',color='b')
plt.xticks(fontsize=f_size)	
plt.yticks(fontsize=f_size)		
plt.ylabel('$CO_2\/(10^{25}\/mol/m^2)$',fontsize=f_size)
ax3 = plt.subplot(313,sharex=ax1)
plt.plot(O2_datetime,O2_col/10**24,linestyle='None',marker='.',color='g')
plt.xticks(fontsize=f_size)	
plt.yticks(fontsize=f_size)		
plt.ylabel('$O_2\/(10^{28}\/mol/m^2)$',fontsize=f_size)

plt.figure(2)
ax5 = plt.subplot(311,sharex=ax1)
plt.plot(O2_datetime,O2_col/10**24,linestyle='None',marker='.',color='g')
plt.xticks(fontsize=f_size)		
plt.yticks(fontsize=f_size)		
plt.ylabel('$O_2\/(10^{28}\/mol/m^2)$',fontsize=f_size)
ax6 = plt.subplot(312,sharex=ax1)
plt.plot(p_datetime,p_pressure,linestyle='-',marker='.',color='k')
plt.xticks(fontsize=f_size)		
plt.yticks(fontsize=f_size)		
plt.ylabel('Pressure (hPa)',fontsize=f_size)
ax7 = plt.subplot(313,sharex=ax1)
plt.plot(O2_datetime,O2_mat['SZA'],linestyle='None',marker='.',color='m')
plt.xticks(fontsize=f_size)
plt.yticks(fontsize=f_size)		
plt.ylabel('SZA (deg)',fontsize=f_size)

plt.show()

