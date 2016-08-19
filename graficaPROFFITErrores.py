import proffitHDF
import matplotlib.pyplot as plt
import numpy as np
import sys

err_dtype = np.dtype([('GAS',np.str_,3), ('BASE_ST',float), ('ILS_ST',float), ('LOS_ST',float), ('SOLAR_ST',float), ('T_ST',float), ('SPEC_ST',float),   ('NOISE_ST',float), ('BASE_SY',float), ('ILS_SY',float), ('LOS_SY',float), ('SOLAR_SY',float), ('T_SY',float), ('SPEC_SY',float), ('NOISE_SY',float)])
err_array = np.empty(2,dtype=err_dtype)


for ii,target_gas in enumerate(['CO2','O2']):
	err_file = '150609.0_072516SN.bin.hdf5'#'150609.44_135143SN.bin.hdf5'#'150801.16_125532SN.bin.hdf5'#'150803.2_075751SN.bin.hdf5'#
	file_path = '/home/D2_RESULTS/PROFFIT_results/ALTZ/'+target_gas+'_vERR/'+err_file.split('.')[0]+'/'
	#file_path = '/home/D2_PROFFIT/Jorge_Files/'+target_gas+'_v1/'+err_file.split('.')[0]+'/'
	err_prf = proffitHDF.proffitHDF(file_path+err_file)

	err_color = ['#d11141', '#00b159', '#00aedb', '#f37735', '#ffc425', '#ffe28a', '#666547']

	print target_gas,file_path,err_file,'SZA =',err_prf.SZA
	print err_prf.especies,err_prf.especies[err_prf.especies == target_gas]
	#print err_prf.VMR[err_prf.especies == target_gas][0],len(err_prf.z)

	VMR_prf = err_prf.VMR[err_prf.especies == target_gas][0]

	BASE_STerr = 100.0*err_prf.BASE_STerr/VMR_prf
	ILS_STerr = 100.0*err_prf.ILS_STerr/VMR_prf
	LOS_STerr = 100.0*err_prf.LOS_STerr/VMR_prf
	SOLAR_STerr = 100.0*err_prf.SOLAR_STerr/VMR_prf
	T_STerr = 100.0*err_prf.T_STerr/VMR_prf
	SPEC_STerr = 100.0*err_prf.SPEC_STerr/VMR_prf
	NOISE_STerr = 100.0*err_prf.NOISE_STerr/VMR_prf

	BASE_SYerr = 100.0*err_prf.BASE_SYerr/VMR_prf
	ILS_SYerr = 100.0*err_prf.ILS_SYerr/VMR_prf
	LOS_SYerr = 100.0*err_prf.LOS_SYerr/VMR_prf
	SOLAR_SYerr = 100.0*err_prf.SOLAR_SYerr/VMR_prf
	T_SYerr = 100.0*err_prf.T_SYerr/VMR_prf
	SPEC_SYerr = 100.0*err_prf.SPEC_SYerr/VMR_prf
	NOISE_SYerr = 100.0*err_prf.NOISE_SYerr/VMR_prf

	print 'PAP ST TOT',np.sqrt(0.249**2 + 0.1749**2 + 0.0249**2 + 0.1049**2 + 0.0449**2 + 0.0249**2)
	print 'PAP SY TOT',np.sqrt(0.17**2 + 0.02**2 + 0.07**2 + 0.01**2 + 0.01**2 + 4.22**2)
	print 'ST',np.mean(BASE_STerr),np.mean(ILS_STerr),np.mean(LOS_STerr),np.mean(SOLAR_STerr),np.mean(T_STerr),np.mean(SPEC_STerr),np.mean(NOISE_STerr)
	print 'ST TOT',np.sqrt(np.mean(BASE_STerr)**2 + np.mean(ILS_STerr)**2 + np.mean(LOS_STerr)**2 + np.mean(SOLAR_STerr)**2 + np.mean(T_STerr)**2 + np.mean(SPEC_STerr)**2 + np.mean(NOISE_STerr)**2)
	print 'SY',np.mean(BASE_SYerr),np.mean(ILS_SYerr),np.mean(LOS_SYerr),np.mean(SOLAR_SYerr),np.mean(T_SYerr),np.mean(SPEC_SYerr)#,np.mean(NOISE_SYerr)
	print 'SY TOT',np.sqrt(np.mean(BASE_SYerr)**2 + np.mean(ILS_SYerr)**2 + np.mean(LOS_SYerr)**2 + np.mean(SOLAR_SYerr)**2 + np.mean(T_SYerr)**2 + np.mean(SPEC_SYerr)**2) #+ np.mean(NOISE_SYerr)**2)

	err_array['GAS'][ii] = target_gas
	err_array['BASE_ST'][ii] = np.mean(BASE_STerr)
	err_array['ILS_ST'][ii] = np.mean(ILS_STerr)
	err_array['LOS_ST'][ii] = np.mean(LOS_STerr)
	err_array['SOLAR_ST'][ii] = np.mean(SOLAR_STerr)
	err_array['T_ST'][ii] = np.mean(T_STerr)
	err_array['SPEC_ST'][ii] = np.mean(SPEC_STerr)
	err_array['NOISE_ST'][ii] = np.mean(NOISE_STerr)
	err_array['BASE_SY'][ii] = np.mean(BASE_SYerr)
	err_array['ILS_SY'][ii] = np.mean(ILS_SYerr)
	err_array['LOS_SY'][ii] = np.mean(LOS_SYerr)
	err_array['SOLAR_SY'][ii] = np.mean(SOLAR_SYerr)
	err_array['T_SY'][ii] = np.mean(T_SYerr)
	err_array['SPEC_SY'][ii] = np.mean(SPEC_SYerr)
	err_array['NOISE_SY'][ii] = np.mean(NOISE_SYerr)

	plt.figure(1)
	plt.suptitle(target_gas+'\nSZA = %.4f\nSTATISTICAL ERRORS' % err_prf.SZA)
	plt.xscale('log')
	plt.plot(BASE_STerr,err_prf.z,label='BASE',color=err_color[0],linewidth=3.0)
	plt.plot(ILS_STerr,err_prf.z,label='ILS',color=err_color[1],linewidth=3.0)
	plt.plot(LOS_STerr,err_prf.z,label='LOS',color=err_color[2],linewidth=3.0)
	plt.plot(SOLAR_STerr,err_prf.z,label='SOLAR',color=err_color[3],linewidth=3.0)
	plt.plot(T_STerr,err_prf.z,label='T',color=err_color[4],linewidth=3.0)
	#plt.plot(SPEC_STerr,err_prf.z,label='SPEC',color=err_color[5],linewidth=3.0) #NO HAY SPEC ERROR EN ST
	plt.plot(NOISE_STerr,err_prf.z,label='NOISE',color=err_color[6],linewidth=3.0)
	#plt.plot(err_prf.VMR[err_prf.especies == target_gas][0],err_prf.z,label='VMR',color='k',linestyle='None',marker='s')
	#plt.plot(err_prf.aPrioriVMR[err_prf.especies == target_gas][0],err_prf.z,label='A PRIORI VMR',color='w',linestyle='None',marker='s')
	plt.legend()
	#plt.savefig('ST_ERR_'+target_gas+'_'+err_file+'.eps')

	plt.figure(2)
	plt.suptitle(target_gas+'\nSZA = %.4f\nSYSTEMATIC ERRORS' % err_prf.SZA)
	plt.xscale('log')
	plt.plot(BASE_SYerr,err_prf.z,label='BASE',color=err_color[0],linewidth=3.0)
	plt.plot(ILS_SYerr,err_prf.z,label='ILS',color=err_color[1],linewidth=3.0)
	plt.plot(LOS_SYerr,err_prf.z,label='LOS',color=err_color[2],linewidth=3.0)
	plt.plot(SOLAR_SYerr,err_prf.z,label='SOLAR',color=err_color[3],linewidth=3.0)
	plt.plot(T_SYerr,err_prf.z,label='T',color=err_color[4],linewidth=3.0)
	plt.plot(SPEC_SYerr,err_prf.z,label='SPEC',color=err_color[5],linewidth=3.0)
	plt.plot(NOISE_SYerr,err_prf.z,label='NOISE',color=err_color[6],linewidth=3.0)
	#plt.plot(err_prf.VMR[err_prf.especies == target_gas][0],err_prf.z,label='VMR',color='k',linestyle='None',marker='s')
	#plt.plot(err_prf.aPrioriVMR[err_prf.especies == target_gas][0],err_prf.z,label='A PRIORI VMR',color='w',linestyle='None',marker='s')
	plt.legend()
	#plt.savefig('SY_ERR_'+target_gas+'_'+err_file+'.eps')

print err_array
#plt.show()
