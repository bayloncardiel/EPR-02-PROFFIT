import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
import datetime as dt

yearAK = 2015
monthAK = 4
dayAK = 26

path_FTIR = ["/home/D2_RESULTS/PROFFIT_results/ALTZ/CO2_v2/%02i%02i%02i/" % (yearAK-2000,monthAK,dayAK),\
"/home/D2_PROFFIT/Jorge_Files/CO2_v1/%02i%02i%02i/" % (yearAK-2000,monthAK,dayAK)]
method_FTIR = ['CO2_v2','CO2_v1']
strategy_FTIR = ['BLOCK','SCALING']

for rr,rlr in enumerate(path_FTIR):
	f_name_FTIR = "%02i%02i%02i.*.hdf5" % (yearAK-2000,monthAK,dayAK)

	list_files = np.array(glob.glob(rlr+f_name_FTIR))
	time_files = np.array([ulu.split("/")[-1].split('_')[-1].split('.')[0][0:-2] for ulu in list_files])
	hour_files = np.array([int(ili[0:2]) for ili in time_files])
	idx_time = np.argsort(time_files)
	time_files = time_files[idx_time]
	hour_files = hour_files[idx_time]
	list_files = list_files[idx_time]
	list_files = list_files[hour_files < 9]
	time_files = time_files[hour_files < 9]

	cm = plt.get_cmap('jet')
	fig = plt.figure(rr)
	ax = fig.add_subplot(111)
	ax.set_color_cycle([cm(1.*i/len(list_files)) for i in range(len(list_files))])
	for oo,ele in enumerate(list_files):
		print 'FILE : %s' % ele.split("/")[-1]
		f = h5py.File(ele,"r")
		for ee in f.keys():					### READING OF RELEVANT VALUES FROM HD5
			if ee == 'CO2':
				for ii in f[ee].keys():
					if ii == method_FTIR[rr]:
						for uu in f[ee][ii].keys():
							if uu == 'NDACC_GMES':
								for aa in f[ee][ii][uu].keys():
									if aa == "CO2.COLUMN_ABSORPTION.SOLAR_AVK":
										col_AK_FTIR = f[ee][ii][uu][aa][()]
									elif aa == 'PRESSURE_INDEPENDENT':
										surfPr_FTIR = f[ee][ii][uu][aa][()]
									elif aa == 'ALTITUDE':
										altLvl_FTIR = f[ee][ii][uu][aa][()]
									elif aa == 'ANGLE.SOLAR_ZENITH.ASTRONOMICAL':
										sza_FTIR = f[ee][ii][uu][aa][()]

		#surfPr_FTIR = np.array([ala/100.0 for ala in surfPr_FTIR])
		#plt.plot(col_AK_FTIR,surfPr_FTIR,linestyle='--',color='r',label='FTIR')
		ax.plot(col_AK_FTIR,altLvl_FTIR,linestyle='-',label='%s:%s:%s' % (time_files[oo][0:2],time_files[oo][2:4],time_files[oo][4:6]))
		plt.xlabel('Averaging Kernels')
		#plt.ylabel('Pressure height (hPa)')
		plt.ylabel('Altitude (km)')
		plt.xlim(0.2,1.4)
		plt.ylim(altLvl_FTIR[0],15.0)
		#plt.gca().invert_yaxis()
		#plt.legend(handlelength=2,fontsize=8)#,ncol=2)
	#plt.suptitle(strategy_FTIR[rr])
	#plt.savefig('AK_150426_06-08_%s.eps' % strategy_FTIR[rr])
plt.show()

