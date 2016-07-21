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

function newlimits_1o(){
##We use here a first order polynomial to obtain the new interval where to look
##for ystable and yunstable
##We need at least to have performed two computations.
local y1smsup y1sminf y2smsup y2sminf
local y1umsup y1uminf y2umsup y2uminf
local ysup yinf
local dist tol step

###We first compute the new value of gamma:
ngamm=$(echo "scale=$digits;($gammini*$ndisty- $gammfin*$disty)/($ndisty- $disty)" | bc)
echo "-1 New gamm=$ngamm, iteration $i"

fpp=1 ##Estimate of the (maximum?) second derivative of the real function in the
      ##interval gammini-gammfin

##The error committed by interpolating is proportional to this value (which we call
##step but has nothing to do with any step):
if [ $(echo "scale=$digits; $ngamm>$gammini && $ngamm < $gammfin" | bc) -eq 1 ]; then
	step=$(echo "scale=$digits; ($ngamm-$gammini)*($ngamm-$gammfin)" | bc) 
	step=$(abs $step)
else
	echo "-1 Problem, ngamm not in the proper interval."
	exit
fi

##Notation for the next iteration:
gamm1=$gammini
gamm2=$ngamm
gamm3=$gammfin
disty1=$disty
disty3=$ndisty
index1=1
index3=2
index2=3

tol=$(echo "scale=$digits; $step*$fpp"| bc)
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

y1minsm=$(echo "scale=$digits;$y1sminf+($ngamm-$gammini)*($y2sminf-$y1sminf)/($gammfin-$gammini) " | bc)
y1maxsm=$(echo "scale=$digits;$y1smsup+($ngamm-$gammini)*($y2smsup-$y1smsup)/($gammfin-$gammini) " | bc)

if [ $(echo "scale=$digits; $y1minsm > $y1maxsm" | bc) -eq 1 ]; then
	aux=$y1minsm
	y1minsm=$y1maxsm
	y1maxsm=$aux
fi
tol=0 ##This is to increase the range
dist=$(echo "scale=$digits; $y1maxsm-$y1minsm" | bc)
y1maxsm=$(echo "scaled=$digits; $y1maxsm+$tol*$dist" | bc)
y1minsm=$(echo "scaled=$digits; $y1minsm-$tol*$dist" | bc)

y1minum=$(echo "scale=$digits;$y1uminf+($ngamm-$gammini)*($y2uminf-$y1uminf)/($gammfin-$gammini) " | bc)
y1maxum=$(echo "scale=$digits;$y1umsup+($ngamm-$gammini)*($y2umsup-$y1umsup)/($gammfin-$gammini) " | bc)
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

function newlimits_1o2(){
##This is a variation of newlimits_1o thought to be continued by a also a first
##order secant method. This is thought to avoid the interpolating errors given by
##the second order interpolation when the interval gets very small.
##Now the zero has to be located between gamm1 and gamm2.

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
if [ $(echo "scale=$digits; $ngamm>$gamm1 && $ngamm < $gamm2" | bc) -eq 1 ]; then
	step=$(echo "scale=$digits; ($ngamm-$gamm1)*($ngamm-$gamm2)" | bc) 
	step=$(abs $step)
else
	echo "-1 Problem, ngamm not in the proper interval."
	exit
fi

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


function newlimits_2o(){
##We use here a second order polynomial to obtain the new interval where to look
##for ystable and yunstable
##We need at least to have performed three computations.
local y1smsup y1sminf y2smsup y2sminf y3smsup y3sminf
local y1umsup y1uminf y2umsup y2uminf y3umsup y3uminf
local ysup yinf
local dist tol
local x1 x2 x3
local y1 y2 y3
local a b c
local nindex1 nindex2 nindex3

##This does no longer uses ordered_gammas, we use here indexes instead.

###We first interpolate. (We need that last three values of gamm and disty)
x1=$gamm1
x2=$gamm2
x3=$gamm3
y1=$disty1
y2=$disty2
y3=$disty3
##Taken from interpolation_2n_order.mw
a=$(echo "scale=100;-($y2*$x1- $y2*$x3- $y1*$x2+$y1*$x3- $y3*$x1+$y3*$x2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
b=$(echo "scale=100;(- $y1*$x2^2+$y1*$x3^2+$y2*$x1^2- $y2*$x3^2- $y3*$x1^2+$y3*$x2^2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
c=$(echo "scale=100; ($x2^2*$x3*$y1- $x2*$x3^2*$y1+$x1^2*$x2*$y3- $x1*$x2^2*$y3- $x1^2*$x3*$y2+$x1*$x3^2*$y2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)

#tmp=$(echo "scale=$digits; $a*$gamm3^2+$b*$gamm3+$c" | bc)
#echo $tmp
#exit
#echo "$a $b $c"
#ngamm=$(echo "scale=$digits; $x2-($x3-$x2)/($y3-$y2)*$y2" | bc)
ngamm=$(echo "scale=$digits; (- $b+sqrt($b^2- 4*$a*$c))/(2*$a)"| bc)
#echo $ngamm
#exit

if [ $(echo "scale=$digits; $ngamm>$x1 && $ngamm < $x2" | bc) -eq 1 ]; then
	gamm1=$x1
	gamm2=$ngamm
	gamm3=$x2
	disty1=$disty1
	disty3=$disty2
	nindex1=$index1
	nindex3=$index2
elif [ $(echo "scale=$digits; $ngamm>$x2 && $ngamm < $x3" | bc) -eq 1 ]; then 
	gamm1=$x2
	gamm2=$ngamm
	gamm3=$x3
	disty1=$disty2
	disty3=$disty3
	nindex1=$index2
	nindex3=$index3
else
	ngamm=$(echo "scale=$digits; (- $b- sqrt($b^2- 4*$a*$c))/(2*$a)"| bc)
	if [ $(echo "scale=$digits; $ngamm>$x1 && $ngamm < $x2" | bc) -eq 1 ]; then
		gamm1=$x1
		gamm2=$ngamm
		gamm3=$x2
		disty1=$disty1
		disty3=$disty2
		nindex1=$index1
		nindex3=$index2
	elif [ $(echo "scale=$digits; $ngamm>$x2 && $ngamm < $x3" | bc) -eq 1 ]; then 
		gamm1=$x2
		gamm2=$ngamm
		gamm3=$x3
		disty1=$disty2
		disty3=$disty3
		nindex1=$index2
		nindex3=$index3
	else
		echo "Error when calculating ngamm, not in the proper interval."
		exit
	fi
fi
echo "-1 New gamm=$ngamm, iteration $i"

fppp=1 ##Estimate of the (maximum?) of the third derivative in the interval
##The error committed by interpolating is proportional to this value (which we call
##step but has nothing to do with any step):
step=$(echo "scale=$digits; ($ngamm- $x1)*($ngamm- $x2)*($ngamm- $x3)" | bc) 
step=$(abs $step)

tol=$(echo "scale=$digits; $step*$fppp"| bc)
get_ys process-stable$index1 $tol
y1smsup=$ysup
y1sminf=$yinf
get_ys process-stable$index2 $tol
y2smsup=$ysup
y2sminf=$yinf
get_ys process-stable$index3 $tol
y3smsup=$ysup
y3sminf=$yinf

get_ys process-unstable$index1 $tol
y1umsup=$ysup
y1uminf=$yinf
get_ys process-unstable$index2 $tol
y2umsup=$ysup
y2uminf=$yinf
get_ys process-unstable$index3 $tol
y3umsup=$ysup
y3uminf=$yinf


y1=$y1sminf
y2=$y2sminf
y3=$y3sminf
a=$(echo "scale=$digits;-($y2*$x1- $y2*$x3- $y1*$x2+$y1*$x3-$y3*$x1+$y3*$x2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
b=$(echo "scale=$digits;(- $y1*$x2^2+$y1*$x3^2+$y2*$x1^2- $y2*$x3^2- $y3*$x1^2+$y3*$x2^2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
c=$(echo "scale=$digits; ($x2^2*$x3*$y1- $x2*$x3^2*$y1+$x1^2*$x2*$y3- $x1*$x2^2*$y3- $x1^2*$x3*$y2+$x1*$x3^2*$y2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
y1minsm=$(echo "scale=$digits; $a*$ngamm^2+$b*$ngamm+$c" | bc)

y1=$y1smsup
y2=$y2smsup
y3=$y3smsup
a=$(echo "scale=$digits;-($y2*$x1- $y2*$x3- $y1*$x2+$y1*$x3-$y3*$x1+$y3*$x2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
b=$(echo "scale=$digits;(- $y1*$x2^2+$y1*$x3^2+$y2*$x1^2- $y2*$x3^2- $y3*$x1^2+$y3*$x2^2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
c=$(echo "scale=$digits; ($x2^2*$x3*$y1- $x2*$x3^2*$y1+$x1^2*$x2*$y3- $x1*$x2^2*$y3- $x1^2*$x3*$y2+$x1*$x3^2*$y2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
y1maxsm=$(echo "scale=$digits; $a*$ngamm^2+$b*$ngamm+$c" | bc)
if [ $(echo "scale=$digits; $y1minsm>$y1maxsm"| bc) -eq 1 ]; then
	aux=$y1minsm
	y1minsm=$y1maxsm
	y1maxsm=$aux
fi
tol=0 ##This is to increase the range
dist=$(echo "scale=$digits; $y1maxsm-$y1minsm" | bc)
y1maxsm=$(echo "scaled=$digits; $y1maxsm+$tol*$dist" | bc)
y1minsm=$(echo "scaled=$digits; $y1minsm-$tol*$dist" | bc)


y1=$y1uminf
y2=$y2uminf
y3=$y3uminf
a=$(echo "scale=$digits;-($y2*$x1- $y2*$x3- $y1*$x2+$y1*$x3-$y3*$x1+$y3*$x2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
b=$(echo "scale=$digits;(- $y1*$x2^2+$y1*$x3^2+$y2*$x1^2- $y2*$x3^2- $y3*$x1^2+$y3*$x2^2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
c=$(echo "scale=$digits; ($x2^2*$x3*$y1- $x2*$x3^2*$y1+$x1^2*$x2*$y3- $x1*$x2^2*$y3- $x1^2*$x3*$y2+$x1*$x3^2*$y2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
y1minum=$(echo "scale=$digits; $a*$ngamm^2+$b*$ngamm+$c" | bc)

y1=$y1umsup
y2=$y2umsup
y3=$y3umsup
a=$(echo "scale=$digits;-($y2*$x1- $y2*$x3- $y1*$x2+$y1*$x3-$y3*$x1+$y3*$x2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
b=$(echo "scale=$digits;(- $y1*$x2^2+$y1*$x3^2+$y2*$x1^2- $y2*$x3^2- $y3*$x1^2+$y3*$x2^2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
c=$(echo "scale=$digits; ($x2^2*$x3*$y1- $x2*$x3^2*$y1+$x1^2*$x2*$y3- $x1*$x2^2*$y3- $x1^2*$x3*$y2+$x1*$x3^2*$y2)/(($x1- $x2)*($x2- $x3)*($x1- $x3))" | bc)
y1maxum=$(echo "scale=$digits; $a*$ngamm^2+$b*$ngamm+$c" | bc)

if [ $(echo "scale=$digits; $y1minum > $y1maxum"| bc) -eq 1 ]; then
	aux=$y1minum
	y1minum=$y1maxum
	y1maxum=$aux
fi
dist=$(echo "scale=$digits; $y1maxum-$y1minum" | bc)
y1maxum=$(echo "scaled=$digits; $y1maxum+$tol*$dist" | bc)
y1minum=$(echo "scaled=$digits; $y1minum-$tol*$dist" | bc)


index1=$nindex1
index2=$i
index3=$nindex3

#echo "$y1umsup $y2umsup $y3umsup $y1maxum"
#echo "$y1uminf $y2uminf $y3uminf $y1minum"
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

tau0=1.687
v0=0.5042
s0=0
delta=0.07

tol=1e-28
tol=$(printf '%0.49f' $tol)
#tol=0.00000000000000000001 #1e-20


#Initial i: ith step to be first computed
i=10

if [ $i -gt 3 ];then
	##In some cases we might have to do the following before continuating
	##with the computations:
	##=====From here
	#ngamm=.511103782577221113928491648107348962237883509797251936664917
	#s=$(echo "scale=$digits; $s0 + $ngamm" | bc)
	#U $ngamm
	#nystable=1.286048501886121250414743327827896444054964699574115230093905
	#nyunstable=1.286048501873594429837832185515165533250420642616936585678173
	#disty2=$(echo "scale=$digits; $nystable- $nyunstable" | bc)
	#echo $ngamm $tau0 $v0 $s0 0 $nystable $u $v $s $disty2 >> diff_stable_unstable.dat
	#exit
	##=====Until here======

	###For the first order methods, the zeros has to be located between
	## gamm1 and gamm2
	gamm1=.791513147808598679756598542600690178586412289493075720943174
	gamm2=.791513147965791897833431621433435141595096386585727083154417
	disty1=.000000000000000000000000000001981465107509923722637289265568
	disty2=-.000000000008296529353095554485574122939165150830042166744787
	index1=9
	index2=5
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
	if [ $(echo "scale=$digits; $nisty*$disty1<0" | bc) -eq 1 ];then
		gamm2=$ngamm
		disty2=$ndisty
		index2=$i
	elif [ $(echo "scale=$digits; $ndisty*$disty2<0" | bc) -eq 1 ];then 
		gamm1=$ngamm
		disty1=$ndisty
		index1=$i
	else 
		echo "-1 Error, ngamm not in the proper interval. Exiting"
		exit
	fi
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
	if [ $(echo "scale=$digits; $ndisty*$disty1<0" | bc) -eq 1 ];then
		gamm2=$ngamm
		disty2=$ndisty
		index2=$i
	elif [ $(echo "scale=$digits; $ndisty*$disty2<0" | bc) -eq 1 ];then 
		gamm1=$ngamm
		disty1=$ndisty
		index1=$i
	else 
		echo "-1 Error, ngamm not in the proper interval. Exiting"
		exit
	fi

	andisty=$(abs $ndisty)
	echo "Error: $andisty"
done

echo "Finished"
