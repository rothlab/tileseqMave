#!/bin/bash

datadir=$1
countdir=${datadir}counts/
paramfile=${datadir}parameters.json
if [[ ! -r $paramfile ]]
	then
	echo "Parameter file is not readable!"
	exit 1
fi

for countfile in `ls ${countdir}counts_sample*.csv`
do
	echo $countfile
	#run in background for parallelization
	Rscript inst/scripts/translateCounts.R $countfile $paramfile>/dev/null&
done
