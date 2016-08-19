import glob
import os
import time

def read_FILESCAN(relativespecpath,token):
	f = open(relativespecpath+'FILESCAN.ERG','r')
	f_lines = f.readlines()
	f.close()
	spc_list = []
	for ii,line in enumerate(f_lines):
		if line.strip() == '$':
			date = f_lines[ii+2].strip()
			if date[0:len(token)] == token:
				num_spectra = int(f_lines[ii+3].strip())
				for ele in f_lines[ii+4:ii+4+num_spectra]:
					spc_list.append(relativespecpath+'%s/cal/' % date+ele.strip())
					
	return spc_list

procmax = 8
pp = 0
ppp = 0
pppsegment = 20
relativespecpath = '/home/D2_HR125/EM27/calpy_EM27/'	###	CARPETA DONDE ESTAN LOS ARCHIVOS *.BIN A USAR
resultpath = '/home/D2_HR125/EM27/RESULTS/'			###	CARPETA DONDE SE GUARDAN LS RESULTADOS
method = 'O2_EM27_v1'								###	GAS QUE SE VA A RECUPERAR, BUSCA EN LA CARPETA DE SKELETON
TOKEN = '16'							###	TOKEN DE FECHAS
station = 'CCA'

fol_path = '/home/D2_PROFFIT/PROFFIT/'

spc_files = read_FILESCAN(relativespecpath,TOKEN)

for files in spc_files:

	files_name = files.split('/')[-1]
	files_path = files.replace(files_name,'')
	files_date = files_name.split('_')[0]
	print 'PROCESSING %s' % (files_name)

	fresult = resultpath+method+'/'+files_date+'/'+files_name+'.hdf5'
	if os.path.isfile(fresult): 
		print '%s EXISTS' % fresult
		continue

	fzipresult = resultpath+method+'/'+files_date+'/'+files_name+'.zip'
	if os.path.isfile(fzipresult):
		print '%s EXISTS' % fzipresult
		continue
	
	while True:
		pp = pp + 1
		pp = pp % procmax
		ppp = pp + pppsegment

		statefile = fol_path+"%i_proffit/retrievalstate.txt" % ppp
		state = int(open(statefile,'r').readline().strip())
		print 'STATE OF %i_proffit: %i' % (ppp,state)
		if state == 0:
			print "RETRIVING IN %i_proffit" % ppp
			print >> open(statefile,'w'), '1'
			time.sleep(0.5)
			sh_str = 'python %srunmethode_pp_v2.py %s %s %s %s %i %s\necho 0 > %s'% (fol_path,method, files_path, files_name, resultpath, ppp, station, statefile)
			start_path = '%sstart_%i.sh' % (fol_path,ppp)
			open(start_path,'w').writelines(sh_str)
			os.system('sh %s &' % start_path)
			break
		print 'WAITING'
		time.sleep(1.0)
			

	
	
	
