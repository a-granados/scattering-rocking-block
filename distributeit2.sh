#!/bin/bash

###This script does the same as distributeit.sh but, instead of passing x1min h1x etc to "2dscan_client.sh", it writes in a file 
###all the points that each client has to compute. Then, 2dscan_client2.sh should read them from the file.
datadir=`pwd`
clientfile=clientlist
clients=`cat $clientfile | grep "^[^#]"`
#clients=`cat $clientfile`

#x1min=-0.1
#xmax=
#xpoints=36 ##This has to be a multiple of the total number of processess!
#hx1=$(echo "scale=10;($x1max - $x1min)/($xpoints-1)" | bc)
x1=0

ypoints=2500

y1min=0.99
y1max=1.01
hy1=$(echo "scale=10; ($y1max - $y1min)/($ypoints-1)" | bc)
#hy1=0.00011

x2=0
y2=0.3

t0=0

delta=0.001

numiterates=50
echo "$delta $numiterates">"values2.tmp"

i=1
for client in $clients; do
	###We should check whether the clients are alive or not
	clients2[$i]=$client
	let i=$i+1
done


let numclients=$i/2
#--Number of processes to be launched:
nprocess=0
for i in `seq 1 1 $numclients`; do
	let index=2*$i
	let nprocess=$nprocess+${clients2[$index]}
done


###I will assume that the total number of points in x is a multiple of the total number of processesses.
let niter=$ypoints/$nprocess ###Number of simulations that each process will have to perform
yminaux=$y1min
for i in `seq 1 1 $numclients`; do
	let clientindex=2*$i-1
	let nprocessindex=2*$i
	nprocessclient=${clients2[$nprocessindex]}
	actualclient=${clients2[$clientindex]}
	y1maxaux=$(echo "scale=10; $yminaux + ($nprocessclient*$niter - 1)*$hy1" | bc)
#-------------We write the points in a file------------
	for y1 in `seq $yminaux $hy1 $y1maxaux`; do
#		for y1 in `seq $y1min $hy1 $y1max`;do
		echo $x1 $y1 $x2 $y2 $t0 >> initial_conditions_$actualclient.dat
#		done
	done
#------------And we launch 2dscan_client2.sh in each client--------------
	(
	ssh $actualclient "cd $datadir && scan_yh2.sh $nprocessclient initial_conditions_$actualclient.dat > process_$actualclient"

	)&
	yminaux=$(echo "scale=10; $yminaux+$nprocessclient*$niter*$hy1" | bc)
done
