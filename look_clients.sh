#!/bin/bash

clientfile=clients_sm
if [ -a clients_sm_tmp ]; then
	clientsfiletmp=clients_sm_tmp
else
	clientsfiletmp=clients_look_tmp
	cat $clientfile | grep "^[^#]" > $clientsfiletmp
fi

numclients=$(cat $clientsfiletmp | wc -l)

for i in `seq 1 1 $numclients`; do
	actualclient=$(tail -n +$i $clientsfiletmp | head -n 1 | awk '{print $1}')
	actcliprocsm=$(ssh $actualclient 'ps -C find_sm -o pid --no-headers | wc -l | tr -d " "'</dev/null)
	actcliprocum=$(ssh $actualclient 'ps -C find_um -o pid --no-headers | wc -l | tr -d " "'</dev/null)
	actcliprocresn=$(ssh $actualclient 'ps -C resn -o pid --no-headers | wc -l | tr -d " "'</dev/null)
	echo "Processes in $actualclient: find_sm=$actcliprocsm find_um=$actcliprocum resn=$actcliprocresn"

done
#rm -f $clientsfiletmp
