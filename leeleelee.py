import numpy as np

f = open('PRFERR_SZA=70.txt','r')
lines = f.readlines()
f.close()

for ii,lin in enumerate(lines):
	if lin.strip() == 'CaF2':
		CO2_err_CaF2 = lines[ii+1].strip().split()
		O2_err_CaF2 = lines[ii+2].strip().split()
	if lin.strip() == 'KBr':
		CO2_err_KBr = lines[ii+1].strip().split()
		O2_err_KBr = lines[ii+2].strip().split()

CO2_err_CaF2 = [float(ele) for ele in CO2_err_CaF2]
O2_err_CaF2 = [float(ele) for ele in O2_err_CaF2]
CO2_err_KBr = [float(ele) for ele in CO2_err_KBr]
O2_err_KBr = [float(ele) for ele in O2_err_KBr]
 
CO2_err_diff = [100*(1-(CO2_err_CaF2[ee]/ele)) for ee,ele in enumerate(CO2_err_KBr)]
O2_err_diff = [100*(1-(O2_err_CaF2[ee]/ele)) for ee,ele in enumerate(O2_err_KBr)]

print 'CO2 DIFF',CO2_err_diff
print 'O2 DIFF',O2_err_diff
