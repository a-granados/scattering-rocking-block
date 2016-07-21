#!/bin/bash

###This is script is called by distribute2.sh, and it is run once on every client.
###It consits on the same program as 2dscan_cleint.sh but it gets the points to be simulated
###from the file "initial_conditions_$clientname.dat

###----architecture detection---

a=$(bash -c 'env | grep DISTRIB')
a=${a:8:8}
b=$(bash -c 'uname -i')
clientarch=$a$b
clientsystem=$clientarch/find_het

nprocess=$1 ##Number of processes to be launched at the same time (this is ncores on other scripts)
#niter=$2 ###Number of simulations that each process will have to perform
pointsfile=$2

firstit=1 ###Number of the starting point. This should be used in case of continuating a killed simulation.

clientname=$(hostname)


#-------We get some values of parameters
delta=$(cat  "values2.tmp" | awk '{ print $1 }')

##Number of iterations to perform by the stroboscopic Poincaré map
numiterates=$(cat  "values2.tmp" | awk '{ print $2 }')


totaliter=$(cat $pointsfile | wc -l | tr -d " ") ##We count the number of points in pointsfile
let totaliter=$totaliter-$firsit+1

for i in `seq $firstit 1 $totaliter`; do
		while [ $(ps -C find_het -o pid --no-headers | wc -l | tr -d " ") -ge $nprocess ]; do
			sleep 4s
		done

		x1=$(tail -n +$i $pointsfile | head -n 1 | awk '{print $1}')
		y1=$(tail -n +$i $pointsfile | head -n 1 | awk '{print $2}')
		x2=$(tail -n +$i $pointsfile | head -n 1 | awk '{print $3}')
		y2=$(tail -n +$i $pointsfile | head -n 1 | awk '{print $4}')
		t0=$(tail -n +$i $pointsfile | head -n 1 | awk '{print $5}')
		actualstate=$(echo "scale=3; $i/$totaliter*100" | bc)
		echo "$actualstate% launched"
		(
			poincfile="stroboscopic_$clientname-$i.dat"
			$clientsystem "$x1" "$y1" "$x2" "$y2" "$t0" "$delta" "$poincfile" "$numiterates"
		)&
done

while [ $(ps -C integrate -o pid --no-headers | wc -l | tr -d " ") -ge 1 ]; do
	sleep 4s
done

#for j in `seq 1 1 $i`; do
#	cat "stroboscopic_$clientnumber-$j.dat" >> "stroboscopic_$clientnumber.dat"
#done
echo "Simulation finished"
