#!/bin/bash

if [ ! -d logs ]; then
	mkdir logs
fi


while read dir; do
	samp=`basename $dir`
	sbatch cellranger_count.sbatch $dir $samp
	
done < fqdirs.txt


