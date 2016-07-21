#!/bin/bash


###This script works different as the other distributeit scripts. It assumes that there exists a list of points to be simulated.
## Then, every 5 seconds or so, the number of processes in each
## client in the clientlist is checked. Then, if it is smaller than the value given in the list, a new process is launched.
datadir=`pwd`
clientfile=clientlist2
#pointsfile="points.dat"
pointsfile=$1
numpoints=$(cat $pointsfile | wc -l)

numiterates=50
delta=0.001

if [ -e actualstate.inf ]; then
	actualpoint=$(cat actualstate.inf | awk '{print $1}')
else
	actualpoint=1
	echo "$actualpoint 0%" > actualstate.inf
fi

while [ $actualpoint -le $numpoints ]; do
	numclients=$(cat $clientfile | wc -l)
	for i in `seq 1 1 $numclients`; do
	actualclient=$(tail -n +$i $clientfile | head -n 1 | awk '{print $1}')
	#echo $actualclient
	nprocess=$(tail -n +$i $clientfile | head -n 1 | awk '{print $2}')
	#echo "Num processes $nprocess"
	if [ "$actualclient" = "pcnot" -o "$actualclient" = "pcmicron" ]; then
		clientsystem=Fedora14i386_x64/find_het
	else
		a=$(ssh $actualclient "bash -c 'env | grep DISTRIB'"</dev/null)
		a=${a:8:8}
		b=$(ssh $actualclient "bash -c 'uname -i'"</dev/null)
		clientarch=$a$b
		clientsystem=$clientarch/find_het
	fi
	actcliproc=$(ssh $actualclient 'ps -C find_het -o pid --no-headers | wc -l | tr -d " "'</dev/null)
	#actcliproc=0
	#echo $actcliproc
	while [ $actcliproc -lt $nprocess ]; do
		if [ $actualpoint -le $numpoints ]; then
		x1=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $1}')
		y1=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $2}')
		x2=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $3}')
		y2=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $4}')
		t0=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $5}')
		poincfile="stroboscopic_$actualpoint-$actualclient.dat"		
		(
		ssh $actualclient "cd $datadir && $clientsystem $x1 $y1 $x2 $y2 $t0 $delta $poincfile $numiterates > trajectory_$actualpoint.dat"</dev/null
		)&
		sleep 1s
		actualstate=$(echo "scale=3; $actualpoint/$numpoints*100" | bc)
		echo "$actualpoint $actualstate% launched" > actualstate.inf
		let actualpoint=$actualpoint+1
		actcliproc=$(ssh $actualclient 'ps -C find_het -o pid --no-headers | wc -l | tr -d " "'</dev/null)
		#prova=$(ssh $actualclient 'ps -C find_het -o pid --no-headers | wc -l | tr -d " "')
		#let actcliproc=$actcliproc+1
		#echo $actcliproc
		else
		actcliproc=$nprocess
		fi
	done
	done
	sleep 5s
done
##Now we wait until all the processes are finished
totalprocs=1
while [ $totalprocs != 0 ]; do
	totalprocs=0
	numclients=$(cat $clientfile | wc -l)
	for i in `seq 1 1 $numclients`; do
	actualclient=$(tail -n +$i $clientfile | head -n 1 | awk '{print $1}')
	actcliproc=$(ssh $actualclient 'ps -C find_het -o pid --no-headers | wc -l | tr -d " "'</dev/null)
	let totalprocs=$totalprocs+$actcliproc
	done
	sleep 2s
done
