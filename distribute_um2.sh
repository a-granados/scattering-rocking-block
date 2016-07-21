#!/bin/bash



###This script works as distributeit3.sh.
datadir=`pwd`
clientfile=clients_um_tmp

pointsfile=$1
numpoints=$(cat $pointsfile | wc -l)

numiterates=$2
delta=$3

actualpoint=1
##This is in case eixam was not included in the clients list
echo "yes" > eixam_um_finished
while [ $actualpoint -le $numpoints ]; do
	numclients=$(cat $clientfile | wc -l)
	for i in `seq 1 1 $numclients`; do
	actualclient=$(tail -n +$i $clientfile | head -n 1 | awk '{print $1}')
	nprocess=$(tail -n +$i $clientfile | head -n 1 | awk '{print $2}')
	if [ "$actualclient" = "eixam" ]; then
		rm -f points-eixam-um
		echo "no" > eixam_um_finished
		j=1
		startindex=$actualpoint
		while [ $j -le $nprocess ];do
			tail -n +$actualpoint $pointsfile | head -n 1 >> points-eixam-um
			let j=$j+1
			let actualpoint=$actualpoint+1
		done
		(
		rsync points-eixam-um agranados@eixam.upc.es:/home/agranados/scattering_map/
		ssh agranados@eixam.upc.es "cd scattering_map && ./distribute_um2.sh points-eixam-um $numiterates $delta $startindex > process_distribute_um2"
		)&
	else
		niceness=$(tail -n +$i $clientfile | head -n 1 | awk '{print $3}')

		a=$(ssh $actualclient "bash -c 'env | grep DISTRIB'"</dev/null)
		a=${a:8:8}
		b=$(ssh $actualclient "bash -c 'uname -i'"</dev/null)
		clientarch=$a$b
		clientsystem=$clientarch/find_um
		actcliproc=$(ssh $actualclient 'ps -C find_um -o pid --no-headers | wc -l | tr -d " "'</dev/null)
		while [ $actcliproc -lt $nprocess ]; do
			if [ $actualpoint -le $numpoints ]; then
			x1=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $1}')
			y1=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $2}')
			x2=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $3}')
			y2=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $4}')
			t0=$(tail -n +$actualpoint $pointsfile | head -n 1 | awk '{print $5}')
			poincfile="stroboscopic_$actualpoint-um.dat"		
			(
			ssh $actualclient "cd $datadir && nice +$niceness $clientsystem $x1 $y1 $x2 $y2 $t0 $delta $poincfile $numiterates "</dev/null
				)&
			sleep 0.5s
			let actualpoint=$actualpoint+1
			actcliproc=$(ssh $actualclient 'ps -C find_um -o pid --no-headers | wc -l | tr -d " "'</dev/null)
			else
			###This is to quit the loop
			actcliproc=$nprocess
			fi
		done
	fi
	done
	done
sleep 5s
##Now we wait until all the processes are finished
totalprocs=1
while [ $totalprocs != 0 ]; do
	totalprocs=0
	numclients=$(cat $clientfile | wc -l)
	for i in `seq 1 1 $numclients`; do
	actualclient=$(tail -n +$i $clientfile | head -n 1 | awk '{print $1}')
	if [ "$actualclient" != "eixam" ];then
		actcliproc=$(ssh $actualclient 'ps -C find_um -o pid --no-headers | wc -l | tr -d " "'</dev/null)
		let totalprocs=$totalprocs+$actcliproc
	fi
	done
	sleep 2s
done

##We now wait until eixam finishes
tmp=$(cat eixam_um_finished)
while [ "$tmp" != "yes" ];do
	tmp=$(cat eixam_um_finished)
done
