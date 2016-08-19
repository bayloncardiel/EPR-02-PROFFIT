#!/bin/bash
procmax=8
pp=0
ppp=0
pppsegment=20
relativespecpath='/home/D2_HR125/calpy_HR125/150426/cal/' #'/home/D2_PROFFIT/ALTZ/spec/bin_v3/SN/2015/'	###	CARPETA DONDE ESTAN LOS ARCHIVOS *.BIN A USAR
resultpath='/home/D2_RESULTS/PROFFIT_results/ALTZ/'			###	CARPETA DONDE SE GUARDAN LS RESULTADOS
method='CO2_v3'							###	GAS QUE SE VA A RECUPERAR, BUSCA EN LA CARPETA DE SKELETON
TOKEN='*'							###	TOKEN DE FECHAS
station='ALTZ'

###	GENERA UNA LISTA DE LOS ARCHIVOS A USAR Y LA ESCRIBE EN listatodo.txt
FILES=$relativespecpath$TOKEN
echo $FILES > listatodo.txt

cd /home/D2_PROFFIT/PROFFIT/

for f in $FILES
do
	echo PROCESSING $(basename $f) ...
	fff=$(basename $f)
	
###	REVISA SI EXISTE EL HDF5 DEL FILE
	fresult="$resultpath$method/*/"${fff}.hdf5
	#echo $fresult
	if [ -e ${fresult} ]
	then
		echo $fresult exist
  		continue
  	fi

###	REVISA SI EXISTE EL ZIP DEL FILE
	fzipresult="$resultpath$method/*/"${fff}.zip
	#echo $fzipresult
	if [ -e ${fzipresult} ]
	then
		echo $fresult exist
  		continue
	fi

	while :
	do
		pp=$(($pp+1))
		pp=$(( $pp  % $procmax ))
		ppp=$(($pppsegment+$pp))

		statefile=$ppp"_proffit/retrievalstate.txt"
		state=$(cat $statefile)
		echo STATE OF $ppp"_proffit:" $state
		if [  "$state" -eq 0  ]
		then
			echo "RETRIVING IN" $ppp"_proffit"
			echo 1 > $statefile 	
			sleep 0.5
			echo   sh runmethode_pp_v2.sh $method $relativespecpath $fff $resultpath $ppp $station > start_$ppp.sh
			echo   "echo 0 > "$statefile >> start_$ppp.sh
			#xterm -e "sh start_"$ppp".sh" &
			sh "start_"$ppp".sh" &
			state=''
			echo DONE
			break
		fi
		echo WAITING
		sleep 1.0
	done
done
