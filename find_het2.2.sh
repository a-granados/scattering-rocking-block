#!/bin/bash

##This is an inproved version of find_he2.sh where we compute the new interval
##where to look for the new values of ystable and yunstable using the
##computations done so far. We use for this the info contained in the files
##proces_stable and process_unstable.

##This is the same program as find_het.sh but we vary gamma instead of tau.
## That is, our parametrization in the unperturbed het manifold is
##(sigma(0),\phi_\U(tau_0+gamma;0,v_0),s0+gamma)

digits=30

#=========================================
#==============function abs ===============
#=========================================

function abs(){
local aux

if [ $(echo "scale=$digits; $1 <0 " | bc) -eq 1 ]; then
	aux=$(echo "scale=$digits; - $1" | bc)
	echo $aux
else
	echo $1
fi

}

#=========================================
#==============function U  ===============
#=========================================

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
	u=$(echo "scale=$digits; ($v0-1)*e($taux)/2- ($v0+1)*e(- $taux)/2 +1" | bc -l)
	v=$(echo "scale=$digits; ($v0-1)*e($taux)/2+($v0+1)*e(- $taux)/2" | bc -l)
else
	u=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2+($v0+1)*e($alphap- $taux)/2 -1" | bc -l)
	v=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2-($v0+1)*e($alphap- $taux)/2" | bc -l)
fi

}

#=========================================
#==============function get_ys ===========
#=========================================

function get_ys(){

local i aux
local iter
local nlines
local error tol
nlines=$(cat $1 | wc -l)
#iter=$2 ##ith iteration of the previous calculations to be considered
tol=$2

i=1
aux=$(tail -n $i $1 | head -n 1 |awk '{print $1}')
j=0
if [ $i -eq 0 ]; then
	error=$(tail -n $i $1 | head -n 1 |awk '{print $2}')
else
	error=0
fi
##We assume error and tol allways positive
##The error commited in the previous calculation has to be larger than the tol
##given by the extrapolation.
while  [ $(echo "scale=$digits; $error< $tol"| bc) -eq 1 ] && [ $i -lt $nlines ]; do
	let i=$i+1
	aux=$(tail -n $i $1 | head -n 1 |awk '{print $1}')
	while [ $aux -ne "0" ] && [ $i -lt $nlines ]; do
		let i=$i+1
		aux=$(tail -n $i $1 | head -n 1 |awk '{print $1}')
	done
	if [ $i -lt $nlines ]; then
		error=$(tail -n $i $1 | head -n 1 | awk '{print $2}')
		j=$i ##Index for the last found error 
	fi
done
if [ $i -ge $nlines ];then
	echo "-1 Warning!! End of file $1 reached! Taking first interval as new interval."
	let i=$j
fi
ysup=$(tail -n $i $1 | head -n 1 |awk '{print $4}')
yinf=$(tail -n $i $1 | head -n 1 |awk '{print $3}')
}


#=========================================
#==============function halfperiod =======
#=========================================

function halfperiod(){
alpha=$(echo "scale=$digits; l((1+$1)/(1-$1))" | bc -l )
}

#=========================================
#==============function launchit =========
#=========================================

function launchit(){
#echo "$u $v $s $y1minsm $y1maxsm $delta"
(
./find_sm.sh $u $v $s $y1minsm $y1maxsm $delta > process_stable
)&
#echo "$u $v $s $y1minum $y1maxum $delta"
(
./find_um.sh $u $v $s $y1minum $y1maxum $delta > process_unstable
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

#=========================================
#==============function newlimits_1o =====
#=========================================

function newlimits_1o(){
##We use here a first order polynomial to obtain the new interval where to look
##for ystable and yunstable
##We need at least to have performed three computations.
local y1smsup y1sminf y2smsup y2sminf
local y1umsup y1uminf y2umsup y2uminf
local ysup yinf
local dist tol

fpp=1 ##Estimate of the second derivative
tol=$(echo "scale=$digits; $step^2*2*$fpp"| bc)
let aux=$i-2
get_ys process-stable$aux $tol
y1smsup=$ysup
y1sminf=$yinf
let aux=$i-1
get_ys process-stable$aux $tol
y2smsup=$ysup
y2sminf=$yinf

let aux=$i-2
get_ys process-unstable$aux $tol
y1umsup=$ysup
y1uminf=$yinf
let aux=$i-1
get_ys process-unstable$aux $tol
y2umsup=$ysup
y2uminf=$yinf

y1minsm=$(echo "scale=$digits; -$y1sminf+2*$y2sminf" | bc)
y1maxsm=$(echo "scale=$digits; -$y1smsup+2*$y2smsup" | bc)
if [ $(echo "scale=$digits; $y1minsm > $y1maxsm" | bc) -eq 1 ]; then
	aux=$y1minsm
	y1minsm=$y1maxsm
	y1maxsm=$aux
fi
tol=0 ##This is to increase the range
dist=$(echo "scale=$digits; $y1maxsm-$y1minsm" | bc)
y1maxsm=$(echo "scaled=$digits; $y1maxsm+$tol*$dist" | bc)
y1minsm=$(echo "scaled=$digits; $y1minsm-$tol*$dist" | bc)

y1minum=$(echo "scale=$digits; -$y1uminf+2*$y2uminf" | bc)
y1maxum=$(echo "scale=$digits; -$y1umsup+2*$y2umsup" | bc)
if [ $(echo "scale=$digits; $y1minum > $y1maxum"| bc) -eq 1 ]; then
	aux=$y1minum
	y1minum=$y1maxum
	y1maxum=$aux
fi
dist=$(echo "scale=$digits; $y1maxum-$y1minum" | bc)
y1maxum=$(echo "scaled=$digits; $y1maxum+$tol*$dist" | bc)
y1minum=$(echo "scaled=$digits; $y1minum-$tol*$dist" | bc)

#echo "$y1smsup $y2smsup $y1maxsm"
#echo "$y1sminf $y2sminf $y1minsm"


}


#=========================================
#==============function newlimits_2o =====
#=========================================

function newlimits_2o(){
##We use here a second order polynomial to obtain the new interval where to look
##for ystable and yunstable
##We need at least to have performed three computations.
local y1smsup y1sminf y2smsup y2sminf y3smsup y3sminf
local y1umsup y1uminf y2umsup y2uminf y3umsup y3uminf
local ysup yinf
local dist tol

fppp=1 ##Estimate of the third derivative at x1
tol=$(echo "scale=$digits; $step^3*6*$fppp"| bc)
let aux=$i-3
get_ys process-stable$aux $tol
y1smsup=$ysup
y1sminf=$yinf
let aux=$i-2
get_ys process-stable$aux $tol
y2smsup=$ysup
y2sminf=$yinf
let aux=$i-1
get_ys process-stable$aux $tol
y3smsup=$ysup
y3sminf=$yinf

let aux=$i-3
get_ys process-unstable$aux $tol
y1umsup=$ysup
y1uminf=$yinf
let aux=$i-2
get_ys process-unstable$aux $tol
y2umsup=$ysup
y2uminf=$yinf
let aux=$i-1
get_ys process-unstable$aux $tol
y3umsup=$ysup
y3uminf=$yinf
y1minsm=$(echo "scale=$digits; $y1sminf-3*$y2sminf+3*$y3sminf" | bc)
y1maxsm=$(echo "scale=$digits; $y1smsup-3*$y2smsup+3*$y3smsup" | bc)
if [ $(echo "scale=$digits; $y1minsm>$y1maxsm"| bc) -eq 1 ]; then
	aux=$y1minsm
	y1minsm=$y1maxsm
	y1maxsm=$aux
fi
tol=0 ##This is to increase the range
dist=$(echo "scale=$digits; $y1maxsm-$y1minsm" | bc)
y1maxsm=$(echo "scaled=$digits; $y1maxsm+$tol*$dist" | bc)
y1minsm=$(echo "scaled=$digits; $y1minsm-$tol*$dist" | bc)

y1minum=$(echo "scale=$digits; $y1uminf-3*$y2uminf+3*$y3uminf" | bc)
y1maxum=$(echo "scale=$digits; $y1umsup-3*$y2umsup+3*$y3umsup" | bc)
if [ $(echo "scale=$digits; $y1minum > $y1maxum"| bc) -eq 1 ]; then
	aux=$y1minum
	y1minum=$y1maxum
	y1maxum=$aux
fi
dist=$(echo "scale=$digits; $y1maxum-$y1minum" | bc)
y1maxum=$(echo "scaled=$digits; $y1maxum+$tol*$dist" | bc)
y1minum=$(echo "scaled=$digits; $y1minum-$tol*$dist" | bc)

#echo "$y1umsup $y2umsup $y3umsup $y1maxum"
#echo "$y1uminf $y2uminf $y3uminf $y1minum"
}

#=================================================
#============Here starts the script ==============
#=================================================

tau0=1.687
v0=0.5042
s0=0
delta=0.07
y1minsm=1
y1maxsm=1.4
y1minum=1
y1maxum=1.4

gammini=0.5
gammfin=1.5
ndata=10
step=$(echo "scale=$digits; ($gammfin-$gammini)/($ndata-1)" | bc)

U $gammini ###This is for the modified function U

s=$(echo "scale=$digits; $s0+$gammini" | bc)

i=1
launchit

ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)

ngamm=$gammini
echo $ngamm $v0 $tau0 $s0 $disty >> diff_stable_unstable.dat
mv process_unstable process-unstable1
mv process_stable process-stable1

##Second iteration: we use the same y1min and y1max
i=2
ngamm=$(echo "scale=$digits;$gammini+$step" | bc)
s=$(echo "scale=$digits; $s0+$ngamm" | bc)
U $ngamm
launchit
ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)
echo $ngamm $v0 $tau0 $s0 $disty >> diff_stable_unstable.dat
mv process_unstable process-unstable$i
mv process_stable process-stable$i

##Thirth iteration, we calculate new y1min and y1max using a first order
##interpolation
i=3
ngamm=$(echo "scale=$digits;$ngamm+$step" | bc)
s=$(echo "scale=$digits; $s0+$ngamm" | bc)
U $ngamm
newlimits_1o

launchit
ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)
echo $ngamm $v0 $tau0 $s0 $disty >> diff_stable_unstable.dat
mv process_unstable process-unstable$i
mv process_stable process-stable$i

###Oju que això ha de comecar ab i=4
for i in `seq 4 1 $ndata`; do
	ngamm=$(echo "scale=$digits;$ngamm+$step" | bc)
	s=$(echo "scale=$digits; $s0+$ngamm" | bc)
	U $ngamm
	newlimits_2o
	launchit
	ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
	yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
	disty=$(echo "scale=$digits; $ystable-$yunstable" | bc)
	echo $ngamm $v0 $tau0 $s0 $disty >> diff_stable_unstable.dat
	mv process_unstable process-unstable$i
	mv process_stable process-stable$i

done
