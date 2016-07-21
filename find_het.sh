#!/bin/bash

digits=30

function U() {
v0=$1
tau=$2
u=$(echo "scale=$digits; ($v0-1)*e($tau)/2+(-$v0-1)*e(-$tau)/2 +1" | bc -l)
v=$(echo "scale=$digits; ($v0-1)*e($tau)/2-(-$v0-1)*e(-$tau)/2" | bc -l)
}

function halfperiod(){
alpha=$(echo "scale=$digits; l((1+$1)/(1-$1))" | bc -l )
}

function launchit(){
(
find_sm.sh $u $v $t0 $y1min $y1max $delta > process_stable
)&
(
find_um.sh $u $v $t0 $y1min $y1max $delta > process_unstable
)&

sleep 2s

totalprocs=2
while [ $totalprocs -gt 0 ]; do
	sleep 5s
	runningstable=$(ps -C find_sm.sh -o pid --no-headers | wc -l | tr -d " ")
	runningunstable=$(ps -C find_um.sh -o pid --no-headers | wc -l | tr -d " ")
	let totalprocs=$runningstable+$runningunstable
done

}

v0=0.3
t0=0
delta=0.01
y1min=0.9
y1max=1.1

halfperiod $v0
step=$(echo "scale=$digits; $alpha/10" | bc)

#tauini=$(echo "scale=$digits; 3*$step" | bc)
tauini=0
U $v0 $tauini


launchit

ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)

ndisty=$disty
ntau=$tauini
echo $v0 $ntau $t0 $disty >> diff_stable_unstable.dat


#while [ $(echo "scale=$digits; $disty*$ndisty >0 "| bc) -eq 1 -a $(echo "scale=digits; $ntau<$alpha" | bc) -eq 1 ]; do
for i in `seq 1 1 10`; do
	tauini=$ntau
	ntau=$(echo "scale=$digits;$tauini+$step" | bc)
	ndisty=$disty
	U $v0 $ntau
	launchit
	ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
	yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
	disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)
	echo $v0 $ntau $t0 $disty >> diff_stable_unstable.dat
done
