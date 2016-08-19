import matplotlib.pyplot as plt
import numpy as np
import math
import datetime as dt

### DESVIACION ESTANDAR Y EN PORCENTAJE DE MEDICIONES CONTIGUAS DE XCO2, CO2 Y O2, PARA LAS COLUMNAS DE CO2 Y O2 SE USAN LOS VALORES FILTRADOS UTILIZADOS PARA CALCULAR XCO2

txtfile = "1212-1512_XCO2.txt"
dt_txt = np.dtype([('fname',np.str_,30),('epochTime',int),('SZA',float),('XCO2',float),('totColCO2',float),('totColO2',float)])
mat_esp = np.genfromtxt(txtfile,dtype = dt_txt,delimiter = " ")
fechas = np.array([(dt.datetime.fromtimestamp(ele)).date() for ele in mat_esp["epochTime"]] )
horas = np.array([(dt.datetime.fromtimestamp(ele)).time() for ele in mat_esp["epochTime"]] )
varias = np.unique(fechas)

string_XCO2 ='XCO2\n \n'
string_CO2 = 'CO2\n \n'
string_O2 = 'O2\n \n'
suma_XCO2 = []
suma_CO2 = []
suma_O2 = []
hrs = []
for unafecha in varias:
	string_XCO2 = string_XCO2 + 'FECHA %02i%02i%02i' % (unafecha.year - 2000,unafecha.month,unafecha.day)+' TOTAL DE MEDICIONES: %d'%len(horas[fechas == unafecha]) +'\n'
	string_CO2 = string_CO2 + 'FECHA %02i%02i%02i' % (unafecha.year - 2000,unafecha.month,unafecha.day)+' TOTAL DE MEDICIONES: %d'%len(horas[fechas == unafecha]) +'\n'
	string_O2 = string_O2 + 'FECHA %02i%02i%02i' % (unafecha.year - 2000,unafecha.month,unafecha.day)+' TOTAL DE MEDICIONES: %d'%len(horas[fechas == unafecha]) +'\n'
	fh = horas[fechas == unafecha]
	ch_XCO2 = mat_esp['XCO2'][fechas == unafecha]
	ch_CO2 = mat_esp['totColCO2'][fechas == unafecha]
	ch_O2 = mat_esp['totColO2'][fechas == unafecha]
	ch_XCO2 = ch_XCO2[fh.argsort()]			### CUIDADO CON EL ORDEN argsort(), ORDENAR EL QUE VA EN FUNCION DEL PRIMERO Y LUEGO EL PRIMERO !!!
	ch_CO2 = ch_CO2[fh.argsort()]			### CUIDADO CON EL ORDEN argsort(), ORDENAR EL QUE VA EN FUNCION DEL PRIMERO Y LUEGO EL PRIMERO !!!
	ch_O2 = ch_O2[fh.argsort()]			### CUIDADO CON EL ORDEN argsort(), ORDENAR EL QUE VA EN FUNCION DEL PRIMERO Y LUEGO EL PRIMERO !!!
	fh = fh[fh.argsort()]
	datehour = np.array([dt.datetime.combine(unafecha,fff) for fff in fh])
	cont_6 = 0
	for y,ele in enumerate(datehour):		
		if y != 0:
			if (((ele - datehour[y-1]).total_seconds() > 180.0) | (cont_6 > 5)):
				string_XCO2 = string_XCO2 +'HORAS USADAS ('+str(len(hrs))+'): '
				string_CO2 = string_CO2 +'HORAS USADAS ('+str(len(hrs))+'): '
				string_O2 = string_O2 +'HORAS USADAS ('+str(len(hrs))+'): '
				for hora in hrs:
					string_XCO2 = string_XCO2 + str(hora.time())+' '
					string_CO2 = string_CO2 + str(hora.time())+' '
					string_O2 = string_O2 + str(hora.time())+' '
				string_XCO2 = string_XCO2 + '\n' + 'MEAN %.4f' % np.mean(suma_XCO2)+'	STD %.4f' % np.std(suma_XCO2) +'	PERCENT %.4f' % (100.0*(np.std(suma_XCO2)/np.mean(suma_XCO2)))+'\n'
				string_CO2 = string_CO2 + '\n' + 'MEAN %.4e' % np.mean(suma_CO2)+'	STD %.4e' % np.std(suma_CO2) +'	PERCENT %.4f' % (100.0*(np.std(suma_CO2)/np.mean(suma_CO2)))+'\n'
				string_O2 = string_O2 + '\n' + 'MEAN %.4e' % np.mean(suma_O2)+'	STD %.4e' % np.std(suma_O2) +'	PERCENT %.4f' % (100.0*(np.std(suma_O2)/np.mean(suma_O2)))+'\n'
				suma_XCO2 = []
				suma_CO2 = []
				suma_O2 = []
				hrs = []
				cont_6 = 0
		suma_XCO2.append(ch_XCO2[y])
		suma_CO2.append(ch_CO2[y])
		suma_O2.append(ch_O2[y])
		hrs.append(ele)
		cont_6 = cont_6 + 1
	
	string_XCO2 = string_XCO2 +'HORAS USADAS ('+str(len(hrs))+'): '
	string_CO2 = string_CO2 +'HORAS USADAS ('+str(len(hrs))+'): '
	string_O2 = string_O2 +'HORAS USADAS ('+str(len(hrs))+'): '
	for hora in hrs:
		string_XCO2 = string_XCO2 + str(hora.time())+' '	
		string_CO2 = string_CO2 + str(hora.time())+' '	
		string_O2 = string_O2 + str(hora.time())+' '	
	string_XCO2 = string_XCO2 + '\n' + 'MEAN %.4f' % np.mean(suma_XCO2)+'	STD %.4f' % np.std(suma_XCO2) +'	PERCENT %.4f' % (100.0*(np.std(suma_XCO2)/np.mean(suma_XCO2)))+'\n\n'
	string_CO2 = string_CO2 + '\n' + 'MEAN %.4e' % np.mean(suma_CO2)+'	STD %.4e' % np.std(suma_CO2) +'	PERCENT %.4f' % (100.0*(np.std(suma_CO2)/np.mean(suma_CO2)))+'\n\n'
	string_O2 = string_O2 + '\n' + 'MEAN %.4e' % np.mean(suma_O2)+'	STD %.4e' % np.std(suma_O2) +'	PERCENT %.4f' % (100.0*(np.std(suma_O2)/np.mean(suma_O2)))+'\n\n'
	suma_XCO2 = []
	suma_CO2 = []
	suma_O2 = []
	hrs = []
	cont_6 = 0
	
fileout = open('DESV_XCO2II.txt','w')	### ARCHIVO TXT CON PROMEDIOS Y DESVIACIONES ESTANDAR DE MEDICIONES CONTIGUAS
fileout.write(string_XCO2)
fileout.close()
fileout = open('DESV_CO2II.txt','w')	### ARCHIVO TXT CON PROMEDIOS Y DESVIACIONES ESTANDAR DE MEDICIONES CONTIGUAS
fileout.write(string_CO2)
fileout.close()
fileout = open('DESV_O2II.txt','w')	### ARCHIVO TXT CON PROMEDIOS Y DESVIACIONES ESTANDAR DE MEDICIONES CONTIGUAS
fileout.write(string_O2)
fileout.close()
