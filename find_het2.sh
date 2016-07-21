#!/bin/bash


##This is the same program as find_het.sh but we vary gamma instead of tau.
## That is, our parametrization in the unperturbed het manifold is
##(sigma(0),\phi_\U(tau_0+gamma;0,v_0),s0+gamma)

digits=30

#function U() {
#gamm=$1
#v0=$2
#tau0=$3
#u=$(echo "scale=$digits; ($v0-1)*e($tau0+$gamm)/2+(-$v0-1)*e(-$tau0-$gamm)/2 +1" | bc -l)
#v=$(echo "scale=$digits; ($v0-1)*e($tau0+$gamm)/2-(-$v0-1)*e(-$tau0-$gamm)/2" | bc -l)
#}

function U() {
###This is a modified function U taken from find_het4.sh. This allows gamm to be greater that alphap.
local gamm
local alphap
local alpha
local i
local taux
local A
gamm=$1
alphap=$(echo "scale=$digits; l((1 + $v0)/(1- $v0))" | bc -l )
alpha=$(echo "scale=$digits; 2*$alphap" | bc )
A=$(echo "scale=$digits; $tau0 + $gamm" | bc)

##Modulus alpha
if [ $(echo "scale=$digits; $A>=0 " | bc) -eq 1 ]; then
	i=1
		while [ $(echo "scale=$digits; $i*$alpha <$A" | bc) -eq 1 ];do
		let i=$i+1
		done
	taux=$(echo "scale=$digits; $A-($i-1)*$alpha" | bc)
else
	i=-1
		while [ $(echo "scale=$digits; $i*$alpha >$A" | bc) -eq 1 ];do
		let i=$i-1
		done
	taux=$(echo "scale=$digits; $A - $i*$alpha" | bc)
fi

	#if [ $(echo "scale=$digits; $A<$alphap " | bc) -eq 1 ]; then
	if [ $(echo "scale=$digits; $taux<$alphap " | bc) -eq 1 ]; then
	u=$(echo "scale=$digits; ($v0-1)*e($taux)/2+(- $v0-1)*e(- $taux)/2 +1" | bc -l)
	v=$(echo "scale=$digits; ($v0-1)*e($taux)/2-(- $v0-1)*e(- $taux)/2" | bc -l)
else
	u=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2-( - $v0-1)*e($alphap- $taux)/2 -1" | bc -l)
	v=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2+(- $v0-1)*e($alphap- $taux)/2" | bc -l)
fi

}



function halfperiod(){
alpha=$(echo "scale=$digits; l((1+$1)/(1-$1))" | bc -l )
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

tau0=2.198
v0=0.856
s0=0
delta=0.05
y1min=0.8
y1max=1.2

#halfperiod $v0
#step=$(echo "scale=$digits; $alpha/10" | bc)

#tauini=$(echo "scale=$digits; 3*$step" | bc)
gammini=3
gammfin=4
ndata=20
step=$(echo "scale=$digits; ($gammfin-$gammini)/$ndata" | bc)

#U $gammini $v0 $tau0
U $gammini ###This is for the modified function U

s=$(echo "scale=$digits; $s0+$gammini" | bc)

launchit

ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)

ndisty=$disty
ngamm=$gammini
echo $ngamm $v0 $tau0 $s0 $disty >> diff_stable_unstable.dat
mv process_unstable process-unstable1
mv process_stable process-stable1

for i in `seq 1 1 $ndata`; do
	gammini=$ngamm
	ngamm=$(echo "scale=$digits;$gammini+$step" | bc)
	ndisty=$disty
	s=$(echo "scale=$digits; $s0+$ngamm" | bc)
	##U $ngamm $v0 $tau0
	U $ngamm
	launchit
	ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
	yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
	disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)
	echo $ngamm $v0 $tau0 $s0 $disty >> diff_stable_unstable.dat
	mv process_unstable process-unstable$i
	mv process_stable process-stable$i

done
