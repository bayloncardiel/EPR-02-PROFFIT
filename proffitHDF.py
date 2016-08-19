import h5py
import numpy as np
import matplotlib.pyplot as plt
from datetime import date,time,datetime
import scipy.constants as sc
import os

class proffitHDF:
	def __init__(self,f_path):
		f_hdf = h5py.File(f_path,"r")
		for ee in f_hdf.keys():
			try:
				f_hdf[ee].dtype
			except:
				self.test = f_hdf[ee].keys()[0]
				self.fecha = datetime.strptime(f_hdf[ee][self.test]["retrieval"]["date"],'%Y-%m-%d').date()
				self.hora = datetime.strptime(f_hdf[ee][self.test]["retrieval"]["time"],'%H:%M:%S').time()
				####################
				###	ESPECIES Y DOFS
				self.especies = np.array(f_hdf[ee][self.test]["retrieval"]["gas"])
				self.DOFS = np.array([ele for ele in f_hdf[ee][self.test]["retrieval"]["dofs"]])		### ARREGLO CON DOFS DE CADA GAS
				####################
				### 	P, T, Z PROFILES		
				self.z = np.array(f_hdf[ee][self.test]["result"]["z"]/1000.0)
				self.P = np.array(f_hdf[ee][self.test]["retrieval"]["P"])
				self.T = np.array(f_hdf[ee][self.test]["retrieval"]["T"])
				self.niveles = len(self.z)
				####################
				### 	VENTANAS ESPECTRALES		
				ventanas = []
				numsdeonda = []
				espectrosMED = []
				espectrosSIM = []
				residuos = []
				for ele in f_hdf[ee][self.test]["spec"].keys():
					if ele[0:5] == "micro":
						ventanas.append(ele)
				ventanas.sort()
				for ele in ventanas:
					numsdeonda.append(f_hdf[ee][self.test]["spec"][ele]["w"])
					espectrosMED.append(f_hdf[ee][self.test]["spec"][ele]["obs"])
					espectrosSIM.append(f_hdf[ee][self.test]["spec"][ele]["sim"])
					residuos.append(f_hdf[ee][self.test]["spec"][ele]["dif"])
				self.ventanas = np.array(ventanas)
				self.numsDeOnda = np.array(numsdeonda)
				self.espectrosMED = np.array(espectrosMED)
				self.espectrosSIM = np.array(espectrosSIM)
				self.residuos = np.array(residuos)
				self.RMSVent = np.array(f_hdf[ee][self.test]["retrieval"]["rms"])
				try:
					self.SNRVent = np.array([np.mean(self.espectrosMED[ii])/ele for ii,ele in enumerate(self.RMSVent)])
				except:
					self.SNRVent = np.array(np.mean(self.espectrosMED)/self.RMSVent)
				####################
				### 	AVERAGING KERNELS
				self.airmass = np.array(f_hdf[ee][self.test]["retrieval"]["airmass"])			### AIRMASS FACTOR CALCULADO COMO RETRIEVAL*.PY
				U = np.zeros([self.niveles,self.niveles],dtype=float)
				Uinv = np.zeros([self.niveles,self.niveles],dtype=float)
				for ii,ele in enumerate(self.airmass):
					U[ii,ii] = ele
					Uinv[ii,ii] = 1.0/ele
				self.AKMat = np.array([ele for ele in f_hdf[ee][self.test]["retrieval"]["AKvmr"]])	### ARREGLO CON MATRICES DE AVERAGING KERNEL RAW PARA CADA GAS
				self.AKPColMat = []
				self.AKTot = []
				for ele in self.AKMat:
					AUinv = np.dot(ele,Uinv) 
					UAUinv = np.dot(U,AUinv)
					self.AKPColMat.append(UAUinv)
					AK_tot = np.dot(UAUinv.transpose(),np.ones([self.niveles],dtype=float))
					self.AKTot.append(AK_tot)
				self.AKPColMat = np.array(self.AKPColMat)			### ARREGLO CON MATRICES DE AVERAGING KERNEL EN UNIDADES DE PARTIAL COLUMN PARA CADA GAS
				self.AKTot = np.array(self.AKTot)				### ARREGLO CON TOTAL AVERAGING KERNEL PARA CADA GAS PARA CADA GAS
				####################
				### 	PERFILES VMR
				self.VMR = np.array([ele for ele in f_hdf[ee][self.test]["retrieval"]["retvmr"]])		### ARREGLO CON VMR DE CADA GAS
				self.aPrioriVMR = np.array([ele for ele in f_hdf[ee][self.test]["retrieval"]["aprvmr"]])	### ARREGLO CON VMR A PRIORI DE CADA GAS
				factor_ND = [(sc.N_A*ele)/((100**3)*(sc.R)*self.T[ii])for ii,ele in enumerate(self.P)]
				self.numDen = []		
				for vol in self.VMR:
					self.numDen.append([ele*vol[ii] for ii,ele in enumerate(factor_ND)])
				self.numDen = np.array(self.numDen)				###	ARREGLO CON NUMBER DENSITY DE CADA GAS		
				####################
				### 	COLUMNAS PARCIALES
				self.partCol = np.array(f_hdf[ee][self.test]["retrieval"]["pcol"])	### PARTIAL COLUMN PARA TARGET GAS
				self.aPrioriCol = np.array(f_hdf[ee][self.test]["retrieval"]["aprpcol"])		### A PRIORI PARTIAL COLUMN PARA TARGET GAS
				self.zUp = np.array(self.z[1:])		
				self.zDown = np.array(self.z[:-1])
				self.totCol = np.array([ele for ele in f_hdf[ee][self.test]["retrieval"]["cols"]])	### ARREGLO DE TOTAL COLUMN PARA CADA GAS
				####################	
				### 	PARAMETROS DE INPUT
				self.latitud = f_hdf[ee][self.test]["NDACC_GMES"]["LATITUDE.INSTRUMENT"]
				self.longitud = f_hdf[ee][self.test]["NDACC_GMES"]["LONGITUDE.INSTRUMENT"]
				self.altitud = f_hdf[ee][self.test]["NDACC_GMES"]["ALTITUDE.INSTRUMENT"]
				self.semiFOV = f_hdf[ee][self.test]["retrieval"]["semifov"]
				self.JD = f_hdf[ee][self.test]["retrieval"]["jd_utm6"]
				self.SZA = f_hdf[ee][self.test]["retrieval"]["sza"]
				self.AZ = f_hdf[ee][self.test]["retrieval"]["az"]
				self.archivoUsado = f_hdf[ee][self.test]["retrieval"]["filename"]
				self.OPD = f_hdf[ee][self.test]["retrieval"]["opd"]
				apddic = ['BOXCAR','TRIAG','HAMMING','BH 3T','BH 4T','NB WEAK','NB MEDIUM','NB STRONG']
				self.apd = apddic[int(f_hdf[ee][self.test]["retrieval"]["apd"])-1]
				####################	
				### 	PARAMETROS DE CALIDAD
				self.epochTime = f_hdf[ee][self.test]["retrieval"]["tepoch"]
				self.meanWShift = f_hdf[ee][self.test]["retrieval"]["meanwshift"]
				self.meanSignal = f_hdf[ee][self.test]["retrieval"]["meansignal"]
				self.meanRMS = f_hdf[ee][self.test]["retrieval"]["meanrms"]
				self.relativeRMS = self.meanRMS/self.meanSignal
				self.targetGas = f_hdf[ee][self.test]["retrieval"]["targetgas"]
				####################	
				### 	ERROR ANALYSIS
				self.npcols_err = f_hdf[ee][self.test]["error"]["errinp"]["npcols"]
				self.pcolTOT_err = f_hdf[ee][self.test]["error"]["partialcolumns"]["err_TOT"]
				self.BASE_STerr = f_hdf[ee][self.test]["error"]["GRP_ST_"+self.targetGas]["BASE"]
				self.ILS_STerr = f_hdf[ee][self.test]["error"]["GRP_ST_"+self.targetGas]["ILS"]
				self.LOS_STerr = f_hdf[ee][self.test]["error"]["GRP_ST_"+self.targetGas]["LOS"]
				self.SOLAR_STerr = f_hdf[ee][self.test]["error"]["GRP_ST_"+self.targetGas]["SOLAR"]
				self.T_STerr = f_hdf[ee][self.test]["error"]["GRP_ST_"+self.targetGas]["T"]
				self.SPEC_STerr = f_hdf[ee][self.test]["error"]["GRP_ST_"+self.targetGas]["SPEC"]
				self.NOISE_STerr = f_hdf[ee][self.test]["error"]["GRP_ST_"+self.targetGas]["NOISE"]
				self.BASE_SYerr = f_hdf[ee][self.test]["error"]["GRP_SY_"+self.targetGas]["BASE"]
				self.ILS_SYerr = f_hdf[ee][self.test]["error"]["GRP_SY_"+self.targetGas]["ILS"]
				self.LOS_SYerr = f_hdf[ee][self.test]["error"]["GRP_SY_"+self.targetGas]["LOS"]
				self.SOLAR_SYerr = f_hdf[ee][self.test]["error"]["GRP_SY_"+self.targetGas]["SOLAR"]
				self.T_SYerr = f_hdf[ee][self.test]["error"]["GRP_SY_"+self.targetGas]["T"]
				self.SPEC_SYerr = f_hdf[ee][self.test]["error"]["GRP_SY_"+self.targetGas]["SPEC"]
				self.NOISE_SYerr = f_hdf[ee][self.test]["error"]["GRP_SY_"+self.targetGas]["NOISE"]
		
		f_hdf.close()

###########################################################################################################################
###########################################################################################################################


