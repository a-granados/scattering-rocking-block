#!/bin/bash

datadir=`pwd`
#clientfile=clientlist2
clientfile=$2
clients=`cat $clientfile | grep "^[^#]"`

echo "Killing clients in 2 seconds, press Ctrl+c to cancel..."
sleep 2s

#----Getting clients-----
i=1
for client in $clients; do
	###We should check whether the clients are alive or not
	clients2[$i]=$client
	let i=$i+1
done

#------Killing clients-----
let numclients=$i/3

#--Number of processes to be launched:

for i in `seq 1 1 $numclients`; do
	let clientindex=3*$i-2
	actualclient=${clients2[$clientindex]}
	#echo $actualclient
#	ssh $actualclient "killall 2dscan_client.sh integrate"
#	ssh $actualclient "killall scan_yh2.sh find_het integrate "  
	ssh $actualclient "killall $1 "  
done
