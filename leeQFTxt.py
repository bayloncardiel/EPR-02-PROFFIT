import numpy as np

f_path = '/home/jorge/OPUS/CONVERSION_2/'
f_txt = 'qflagsSN_GOOD_2013.txt'
f_QF = open(f_path+f_txt,'r')
lines_QF = f_QF.readlines()
f_QF.close()

fechas = []
for line in lines_QF:
	line_split = line.strip().split()
	if line_split[8] == 'CaF2':
		fechas.append(line_split[0])
fechas = np.array([ele[0:6] for ele in fechas])
str_fechas = ''
for fec in np.unique(fechas):
	str_fechas = str_fechas+fec+'\n'
f_fechas = open('fechasCaF2.txt','a')
f_fechas.write(str_fechas)
f_fechas.close()
print np.unique(fechas)


