#LEE ARCHIVO BIN ESPECIFICADO Y GRAFICA ESPECTRO

import numpy as np
import matplotlib.pyplot as plt
import sys
import struct
import os
import matplotlib.dates as mdates
import datetime as dt

f_path_1 = sys.argv[1]
f_path_2 = sys.argv[2]

for ll,f_path in enumerate([f_path_1,f_path_2]):
	binfile = open(f_path,'r')
	ft = f_path.split('/')[-1]
	date = ft[0:6]
	hour = ft.split('_')[1][0:6]

	print 'DATE & HOUR',dt.datetime.strptime(str(date+'-'+hour),'%y%m%d-%H%M%S')

	contador = 0				###   CONTADOR DE $ EN ARCHIVO PARA CARACTERES
	contador2 = 0				###   CONTADOR DE $ EN ARCHIVO PARA LINEAS
	binstring = ''
	binline = []

	for line in binfile:
		binline.append(line)
		binstring = binstring + line
	binfile.close()

	for k,line in enumerate(binline):
		if line.strip() == '$':
			contador2 = contador2 + 1
			if contador2 == 1:
				print 'ELEVATION', np.pi*np.abs(90.0-float(binline[k+4]))/180.0
				print 'AZIMUTHAL', np.pi*float(binline[k+5])/180.0

	for j,character in enumerate(binstring):
		if character == '$':
			contador = contador + 1
			if contador == 6:
				binnpt = binstring[j+3:j+23]			###   20 BITS DESPUES DE $\N
				binspec = binstring[j+23:len(binstring)]	###   TODOS LOS BITS DESPUES DE BINNPT

	paramsnumonda = struct.unpack('2d1i',binnpt)				###   DOS DOBLES (FIRST_NUE, DELTA_NUE) Y UN ENTERO (NPT)
	spec = struct.unpack(str(paramsnumonda[2])+'f',binspec)			###   NPT FLOATS (ESPECTRO)
	w=np.arange(paramsnumonda[2])*paramsnumonda[1]+paramsnumonda[0]

	if ll == 0:
		ax1 = plt.subplot('21%i' % (ll+1))
	else:	
		plt.subplot('21%i' % (ll+1),sharex=ax1)		
	plt.title(f_path.split('/')[-4])
	plt.plot(w,spec,'b',label=dt.datetime.strptime(str(date+'-'+hour),'%y%m%d-%H%M%S'))
	plt.xlim([6100,6500.0])
	plt.ylim([min(spec),max(spec)])
	plt.legend()

plt.show()

