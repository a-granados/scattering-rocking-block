#!/bin/bash


##This is to be launched on a single computer. Use scan_yh2.sh and distributeit2.sh to distribute the computations
##on several computers.

####!!!!!!!!!!!!!
##This no longer works with the new find_het and integrate. The inputs have to adapted.

nprocess=30 ##Number of processes to be launched at the same time



yhmin=0.9
yhmax=1.1
hyh=0.00011

#numiterates=1 ##Number of iterations to perform by the stroboscopic Poincaré map
i=0
for yh in `seq $yhmin $hyh $yhmax`; do
		sleep 1s
		while [ $(ps -C find_het -o pid --no-headers | wc -l | tr -d " ") -ge $nprocess ]; do
			sleep 0.5s
		done
		let i=$i+1
		(
			outfile="outfile_$i.dat"
			find_het "$yh" "$outfile"

		)&
done

while [ $(ps -C find_het -o pid --no-headers | wc -l | tr -d " ") -ge 1 ]; do
	sleep 0.5s
done


for j in `seq 1 1 $i`; do
	cat "outfile_$j.dat" >> "diffH_yh.dat"
done

rm -f outfile_*

echo "Finished"
