#!/bin/bash


##This program is to continuate a previous execution of find_het4.2.sh

## Our parametrization in the unperturbed het manifold is
##(sigma(0),\phi_\U(tau_0+gamma;0,v_0),s0+gamma)


function U() {
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
#echo $taux

#if [ $(echo "scale=$digits; $A<$alphap " | bc) -eq 1 ]; then
if [ $(echo "scale=$digits; $taux<$alphap " | bc) -eq 1 ]; then
	u=$(echo "scale=$digits; ($v0-1)*e($taux)/2+(- $v0-1)*e(- $taux)/2 +1" | bc -l)
	v=$(echo "scale=$digits; ($v0-1)*e($taux)/2-(- $v0-1)*e(- $taux)/2" | bc -l)
else
	u=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2-( - $v0-1)*e($alphap- $taux)/2 -1" | bc -l)
	v=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2+(- $v0-1)*e($alphap- $taux)/2" | bc -l)
fi

}

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
##The error commited in the previous calculation to locate ystable and yunstable has to be larger than the tol
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


function newlimits_1o2(){
##This is a variation of newlimits_1o thought to be continued by a also a first
##order secant method. This is thought to avoid the interpolating errors given by
##the second order interpolation when the interval gets very small.
##Unlike in find_het5.2_linear.sh, we do not use values with different signs to
##find the new ngamm but just the last two ones. This could improve the results
##because, although we are extrapolating, the last two computed values are closer
##to the real zero.

##We use here a first order polynomial to obtain the new interval where to look
##for ystable and yunstable
##We need at least to have performed two computations.
local y1smsup y1sminf y2smsup y2sminf
local y1umsup y1uminf y2umsup y2uminf
local ysup yinf
local dist tol step


###We first compute the new value of gamma:
ngamm=$(echo "scale=$digits;($gamm1*$disty2- $gamm2*$disty1)/($disty2- $disty1)" | bc)
echo "-1 New gamm=$ngamm, iteration $i"

fpp=1 ##Estimate of the (maximum?) second derivative of the real function in the
      ##interval gammini-gammfin

##The error committed by interpolating is proportional to this value (which we call
##step but has nothing to do with any step):
step=$(echo "scale=$digits; ($ngamm-$gamm1)*($ngamm-$gamm2)" | bc) 
step=$(abs $step)

tol=$(echo "scale=$digits; $step*$fpp"| bc)
let aux=$index1
get_ys process-stable$aux $tol
y1smsup=$ysup
y1sminf=$yinf
let aux=$index2
get_ys process-stable$aux $tol
y2smsup=$ysup
y2sminf=$yinf

let aux=$index1
get_ys process-unstable$aux $tol
y1umsup=$ysup
y1uminf=$yinf
let aux=$index2
get_ys process-unstable$aux $tol
y2umsup=$ysup
y2uminf=$yinf

y1minsm=$(echo "scale=$digits;$y1sminf+($ngamm-$gamm1)*($y2sminf-$y1sminf)/($gamm2-$gamm1) " | bc)
y1maxsm=$(echo "scale=$digits;$y1smsup+($ngamm-$gamm1)*($y2smsup-$y1smsup)/($gamm2-$gamm1) " | bc)

if [ $(echo "scale=$digits; $y1minsm > $y1maxsm" | bc) -eq 1 ]; then
	aux=$y1minsm
	y1minsm=$y1maxsm
	y1maxsm=$aux
fi
tol=0 ##This is to increase the range
dist=$(echo "scale=$digits; $y1maxsm-$y1minsm" | bc)
y1maxsm=$(echo "scaled=$digits; $y1maxsm+$tol*$dist" | bc)
y1minsm=$(echo "scaled=$digits; $y1minsm-$tol*$dist" | bc)

y1minum=$(echo "scale=$digits;$y1uminf+($ngamm-$gamm1)*($y2uminf-$y1uminf)/($gamm2-$gamm1) " | bc)
y1maxum=$(echo "scale=$digits;$y1umsup+($ngamm-$gamm1)*($y2umsup-$y1umsup)/($gamm2-$gamm1) " | bc)
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

function launchit(){
(
find_sm.sh $u $v $s $y1minsm $y1maxsm $delta > process_stable
)&
(
find_um.sh $u $v $s $y1minum $y1maxum $delta > process_unstable
)&

sleep 2s
local totalprocs
totalprocs=2
while [ $totalprocs -gt 0 ]; do
	sleep 5s
	runningstable=$(ps -C find_sm.sh -o pid --no-headers | wc -l | tr -d " ")
	runningunstable=$(ps -C find_um.sh -o pid --no-headers | wc -l | tr -d " ")
	let totalprocs=$runningstable+$runningunstable
done

}

function max(){
	if [ $(echo "scale=$digits; $1>$2" | bc) -eq 1 ]; then
	echo $1
	else
	echo $2
	fi
}

function min(){
	if [ $(echo "scale=$digits; $1>$2" | bc) -eq 1 ]; then
	echo $2
	else
	echo $1
	fi
}

function abs(){
local aux

if [ $(echo "scale=$digits; $1 <0 " | bc) -eq 1 ]; then
	aux=$(echo "scale=$digits; - $1" | bc)
	echo $aux
else
	echo $1
fi

}

#========================================================================
#=================And here starts the script=============================
#========================================================================
digits=60

tau0=0.62758
v0=0.48
s0=0
delta=0.01

tol=1e-26
tol=$(printf '%0.49f' $tol)
#tol=0.00000000000000000001 #1e-20


#Initial i: ith step to be first computed
i=5

if [ $i -gt 3 ];then
	##In some cases we might have to do the following before continuating
	##with the computations:
	##=====From here
	#ngamm=2.454807914781370617746678735711364334375686762458130588857298
	#s=$(echo "scale=$digits; $s0 + $ngamm" | bc)
	#U $ngamm
	#nystable=1.124868898141598916270349054818737950018135039734759691412054
	#nyunstable=1.124868890420583489429050270585745173876990351737641745463498
	#disty2=$(echo "scale=$digits; $nystable- $nyunstable" | bc)
	#echo $ngamm $tau0 $v0 $s0 0 $nystable $u $v $s $disty2 >> diff_stable_unstable.dat
	#exit
	##=====Until here======

	###For the first order methods, the zeros has to be located between
	## gamm1 and gamm2
	##Unlike in find_het5.2_linear.sh, now these values will be just the
	##last ones of the previous execution.
	gamm1=2.454460864663766661763159610629896042226790055566927516051367
	gamm2=2.454807914781370617746678735711364334375686762458130588857298
	disty1=-.000007076769435381537088882932267290237224752722806652778734
	disty2=.000000007721015426841298784232992776141144687997117945948556
	let index1=$i-2
	let index2=$i-1
	newlimits_1o2 ##ngamm is computed in this function
	s=$(echo "scale=$digits; $s0 + $ngamm" | bc)
	U $ngamm
	launchit
	nystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
	nyunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
	ndisty=$(echo "scale=$digits; $nystable- $nyunstable" | bc)
	echo $ngamm $tau0 $v0 $s0 0 $nystable $u $v $s $ndisty >> diff_stable_unstable.dat
	mv process_unstable process-unstable$i
	mv process_stable process-stable$i
	gamm1=$gamm2
	gamm2=$ngamm
	disty1=$disty2
	disty2=$ndisty
	index1=$index2
	index2=$i	
fi


andisty=$(abs $ndisty)
echo "Error: $andisty"
while [ $(echo "scale=$digits; $andisty>$tol "| bc) -eq 1 ];do
	let i=i+1
	newlimits_1o2 ##ngamm is computed in this function
	s=$(echo "scale=$digits; $s0 + $ngamm" | bc)
	U $ngamm
	launchit
	nystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
	nyunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
	ndisty=$(echo "scale=$digits; $nystable- $nyunstable" | bc)
	echo $ngamm $tau0 $v0 $s0 0 $nystable $u $v $s $ndisty >> diff_stable_unstable.dat
	mv process_unstable process-unstable$i
	mv process_stable process-stable$i
	gamm1=$gamm2
	gamm2=$ngamm
	disty1=$disty2
	disty2=$ndisty
	index1=$index2
	index2=$i	

	andisty=$(abs $ndisty)
	echo "Error: $andisty"
done

echo "Finished"
