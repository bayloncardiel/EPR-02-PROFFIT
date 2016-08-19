#!/bin/sh

### GENERA ARCHIVOS HDF USANDO UNA RUTINA DETERMINADA (hdffromretrieval*.py)
FECHA=$1
TEST=$2
RESULTSPATH="/home/D2_RESULTS/PROFFIT_results/ALTZ"
TOKEN="*"

for DATIL in $RESULTSPATH/$TEST/$FECHA/$TOKEN
do
	if [ "${DATIL##*.}" = "zip" ]
	then
		cd $RESULTSPATH/$TEST/$FECHA
		unzip $DATIL
		python /home/D2_PROFFIT/PROFFIT/hdffromretrieval3.py ${DATIL%.*} $TEST
		rm -r ${DATIL%.*}
	fi
done
