#!/bin/bash



###This script works as distributeit3.sh.
datadir=`pwd`
clientfile=clients_sm_tmp

pointsfile=$1
numpoints=$(cat $pointsfile | wc -l)

numiterates=$2
delta=$3

actualpoint=1
while [ $actualpoint -le $numpoints ]; do
	numclients=$(cat $clientfile | wc -l)
	for i in `seq 1 1 $numclients`; do
	actualclient=$(tail -n +$i $clientfile | head -n 1 | awk '{print $1}')
	#echo $actualclient
	nprocess=$(tail -n +$i $clientfile | head -n 1 | awk '{print $2}')
	niceness=$(tail -n +$i $clientfile | head -n 1 | awk '{print $3}')
	#echo "Num processes $nprocess"
	clientsystem=./find_sm
	actcliproc=$(ssh $actualclient 'ps -C find_sm -o pid --no-headers | wc -l | tr -d " "'</dev/null)
	#actcliproc=0
	#echo $actcliproc
	while [ $actcliproc -lt $nprocess ]; do
		if [ $actualpoint -le $numpoints ]; then
		x1=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $1}')
		y1=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $2}')
		x2=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $3}')
		y2=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $4}')
		t0=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $5}')
		#echo $x1 $y1 $x2 $y2 $t0
		#poincfile="stroboscopic_$actualpoint-$actualclient.dat"		
		poincfile="stroboscopic_$actualpoint-sm.dat"		
		(
		#ssh $actualclient "pkill -o resn && cd $datadir && nice -n $niceness $clientsystem $x1 $y1 $x2 $y2 $t0 $delta $poincfile $numiterates && nice -n $niceness resn&"</dev/null
		ssh $actualclient "cd $datadir && nice -n $niceness $clientsystem $x1 $y1 $x2 $y2 $t0 $delta $poincfile $numiterates"</dev/null
		)&
		sleep 0.5s
		let actualpoint=$actualpoint+1
		#actualstate=$(echo "scale=3; $actualpoint/$numpoints*100" | bc)
		#echo "$actualpoint $actualstate% launched" > actualstate.inf
		actcliproc=$(ssh $actualclient 'ps -C find_sm -o pid --no-headers | wc -l | tr -d " "'</dev/null)
		#prova=$(ssh $actualclient 'ps -C find_het -o pid --no-headers | wc -l | tr -d " "')
		#let actcliproc=$actcliproc+1
		#echo $actcliproc
		else
		###This is to quit the loop
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
	actcliproc=$(ssh $actualclient 'ps -C find_sm -o pid --no-headers | wc -l | tr -d " "'</dev/null)
	let totalprocs=$totalprocs+$actcliproc
	done
	sleep 2s
done
