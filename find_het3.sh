#!/bin/bash


##This is the same program as find_het2.sh but using a bolzanos method.

## Our parametrization in the unperturbed het manifold is
##(sigma(0),\phi_\U(tau_0+gamma;0,v_0),s0+gamma)

digits=30

function U() {
gamm=$1
v0=$2
tau0=$3
alphap=$(echo "scale=$digits; l((1+$v0)/(1-$v0))" | bc -l )
alpha=$(echo "scale=$digits; 2*$alphap" | bc )
A=$(echo "scale=$digits; $tau0+$gamm" | bc)

##Modulus alpha
if [ $(echo "scale=$digits; $A>=0 " | bc) -eq 1 ]; then
	local i=1
		while [ $(echo "scale=$digits; $i*$alpha <$A" | bc) -eq 1 ];do
		let i=$i+1
		done
	local taux=$(echo "scale=$digits; $A-($i-1)*$alpha" | bc)
else
	local i=-1
		while [ $(echo "scale=$digits; $i*$alpha >$A" | bc) -eq 1 ];do
		let i=$i-1
		done
	local taux=$(echo "scale=$digits; $A-$i*$alpha" | bc)
fi


if [ $(echo "scale=$digits; $A<$alphap " | bc) -eq 1 ]; then
	u=$(echo "scale=$digits; ($v0-1)*e($taux)/2+(-$v0-1)*e(-$taux)/2 +1" | bc -l)
	v=$(echo "scale=$digits; ($v0-1)*e($taux)/2-(-$v0-1)*e(-$taux)/2" | bc -l)
else
	u=$(echo "scale=$digits; -($v0-1)*e(-$alphap+$taux)/2-(-$v0-1)*e($alphap-$taux)/2 -1" | bc -l)
	v=$(echo "scale=$digits; -($v0-1)*e(-$alphap + $taux)/2+(-$v0-1)*e($alphap-$taux)/2" | bc -l)
fi

}


function launchit(){
(
find_sm.sh $u $v $s $y1min $y1max $delta > process_stable
)&
(
find_um.sh $u $v $s $y1min $y1max $delta > process_unstable
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

function newlimits_y1(){
###I assume that ystable and yunstable are both positive
if [ $(echo "scale=$digits; $ystable>$yunstable" | bc) -eq 1 ]; then
	y1max=$(echo "scale=$digits; $ystable+2*($ystable-$yunstable)" | bc)
	y1min=$(echo "scale=$digits; $yunstable-2*($ystable-$yunstable)" | bc)
else
	y1max=$(echo "scale=$digits; $yunstable-2*($ystable-$yunstable)" | bc)
	y1min=$(echo "scale=$digits; $ystable+2*($ystable-$yunstable)" | bc)
fi

}

tau0=0
v0=0.3
s0=0
delta=0.01
y1min=0.9
y1max=1.1


##We assume that the diff function has different 
##signs at these values
gammini=2.34
gammfin=2.42

##Put some values which agree with the sign of dist at gammaini and gammfin:
distini=1
distfin=-1
distgamm=$(echo "scale=$digits; $gammfin-$gammini" | bc)

tol=0.0000000000001


while [ $(echo "scale=$digits; $distgamm>$tol "| bc) -eq 1 ];do
	ngamm=$(echo "scale=$digits; $gammini+($gammfin-$gammini)/2" | bc)
	s=$(echo "scale=$digits; $s0+$ngamm" | bc)
	U $ngamm $v0 $tau0
	launchit
	ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
	yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
	disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)
	echo $ngamm $tau0 $v0 $s0 0 $ystable $u $v $s $disty >> diff_stable_unstable.dat
	if [ $(echo "scale=$digits; $disty*$distini<0" | bc) -eq 1 ]; then
		gammfin=$ngamm
	else
		gammini=$ngamm
	fi
	distgamm=$(echo "scale=$digits; $gammfin-$gammini" | bc)
	echo " $distgamm"
	newlimits_y1 ###This is very dangerous. We try to solve it in find_het4
done
