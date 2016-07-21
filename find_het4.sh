#!/bin/bash


##This is the same program as find_het3.sh but instead of bolzanos method we use here the secant method.
###For the recalulation of the newlimits, ymax and ymin, we take the maximum of the two previous one.

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
echo $taux

#if [ $(echo "scale=$digits; $A<$alphap " | bc) -eq 1 ]; then
if [ $(echo "scale=$digits; $taux<$alphap " | bc) -eq 1 ]; then
	u=$(echo "scale=$digits; ($v0-1)*e($taux)/2+(- $v0-1)*e(- $taux)/2 +1" | bc -l)
	v=$(echo "scale=$digits; ($v0-1)*e($taux)/2-(- $v0-1)*e(- $taux)/2" | bc -l)
else
	u=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2-( - $v0-1)*e($alphap- $taux)/2 -1" | bc -l)
	v=$(echo "scale=$digits; -($v0-1)*e(- $alphap + $taux)/2+(- $v0-1)*e($alphap- $taux)/2" | bc -l)
fi

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

function newlimits_y1(){
###I assume that ystable and yunstable are both positive
local aux1
local aux2
aux1=$(max $ystable $yunstable)
aux2=$(max $nystable $nystable)
y1max=$(max $aux1 $aux2)
aux1=$(min $ystable $yunstable)
aux2=$(min $nystable $nyunstable)
y1min=$(min $aux1 $aux2)
}

function newlimits2_y1(){
##In this version we find in diff_stable_unstable.dat the interval
##where ngam is located. Then we look for y1min and y1max in unstable_manifold.dat
##and stable_manifold.dat
##We compute y1max and y1min for the stable and the unstable manifolds, y1maxsm, y1minsm,...
##We should have calculated at least two iterates.

local numdata
local file1

local k
local aux


file1=ordered_gammas

numdata=$(cat $file1.dat | wc -l)

aux=$(tail -n +1 $file1.dat | head -n 1 | awk '{print $1}')

if [ $(echo "scale=$digits; $ngamm<$aux" | bc) -eq 1 ]; then
	echo "$ngamm $i" > $file1.tmp
	y1maxsm=$firsty1maxsm
	y1minsm=$firsty1minsm
	y1maxum=$firsty1maxum
	y1minum=$firsty1minum
	cat $file1.tmp $file1.dat > $file1.tmp2
	mv -f $file1.tmp2 $file1.dat
	rm -f $file1.tmp

else
	k=1
	while [ $k -le $numdata ] && [ $(echo "scale=$digits; $aux<$ngamm" | bc) -eq 1 ]; do
		let k=k+1
		aux=$(tail -n +$k $file1.dat | head -n 1 | awk '{print $1}')
	done

	echo "$ngamm $i" > $file1.tmp
	if [ $k -gt $numdata ]; then
		y1maxsm=$firsty1maxsm
		y1minsm=$firsty1minsm
		y1maxum=$firsty1maxum
		y1minum=$firsty1minum
		cat $file1.dat $file1.tmp > $file1.tmp2
		mv -f $file1.tmp2 $file1.dat
		rm -f $file1.tmp
	else
		local index
		local y1
		local y2
		local y3
		local y4
		local aux2
		local aux3
		index=$(tail -n +$k $file1.dat | head -n 1 | awk '{print $2}')
		y1=$(tail -n +$index unstable_manifold.dat  | head -n 1 | awk '{print $2}')
		y2=$(tail -n +$index stable_manifold.dat  | head -n 1 | awk '{print $2}')
		let k=k-1
		index=$(tail -n +$k $file1.dat | head -n 1 | awk '{print $2}')
		y3=$(tail -n +$index unstable_manifold.dat  | head -n 1 | awk '{print $2}')
		y4=$(tail -n +$index stable_manifold.dat  | head -n 1 | awk '{print $2}')
#		aux2=$(max $y1 $y2)
#		aux3=$(max $y3 $y4)
#		y1max=$(max $aux2 $aux3)
#		aux2=$(min $y1 $y2)
#		aux3=$(min $y3 $y4)
#		y1min=$(min $aux2 $aux3)
		y1maxsm=$(max $y2 $y4)
		y1minsm=$(min $y2 $y4)
		y1maxum=$(max $y1 $y3)
		y1minum=$(min $y1 $y3)
		head -n $k $file1.dat > $file1.tmp2
		let aux=$numdata-$k
		tail -n $aux $file1.dat > $file1.tmp3
		cat $file1.tmp2 $file1.tmp $file1.tmp3 > $file1.dat
		rm -f $file1.tmp*
	fi
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

########
##############And here starts the script###############
#######
digits=60

tau0=0.42
v0=0.4
s0=0
delta=0.2

##These are the initial values where to look for yu and ys
firsty1minsm=0.9
firsty1maxsm=1.1
firsty1minum=0.9
firsty1maxum=1.1
y1maxsm=$firsty1maxsm
y1minsm=$firsty1minsm
y1maxum=$firsty1maxum
y1minum=$firsty1minum

##We will need that the diff function has different 
##signs at these values. Well, it is no longer needed with the secant method,
##but it would help, a lot.
gammini=0.76
gammfin=0.79
tol=0.00000000000000001
##First iteration:
i=1
echo "$gammini $i" > ordered_gammas.dat
s=$(echo "scale=$digits; $s0 + $gammini" | bc)
U $gammini
launchit
ystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
yunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
mv process_unstable process-unstable$i
mv process_stable process-stable$i
disty=$(echo "scale=$digits; $ystable- $yunstable" | bc)
echo $gammini $tau0 $v0 $s0 0 $ystable $u $v $s $disty >> diff_stable_unstable.dat

##Second iteration
let i=i+1
echo "$gammfin $i" >> ordered_gammas.dat

###Outdated
##This is dangerous:
#aux1=$(max $ystable $yunstable)
#aux2=$(min $ystable $yunstable)
##
#y1max=$(echo "scale=$digits; $aux1+2*($aux1- $aux2)" | bc)
#y1min=$(echo "scale=$digits; $aux2-2*($aux1- $aux2)" | bc)
##

s=$(echo "scale=$digits; $s0+ $gammfin" | bc)
U $gammfin
launchit
nystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
nyunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
ndisty=$(echo "scale=$digits; $nystable- $nyunstable" | bc)
echo $gammfin $tau0 $v0 $s0 0 $nystable $u $v $s $ndisty >> diff_stable_unstable.dat

mv process_unstable process-unstable$i
mv process_stable process-stable$i

andisty=$(abs $ndisty)
echo "Error: $andisty"
while [ $(echo "scale=$digits; $andisty>$tol "| bc) -eq 1 ];do
	let i=i+1
	ngamm=$(echo "scale=$digits;($gammini*$ndisty- $gammfin*$disty)/($ndisty- $disty)" | bc)
	newlimits2_y1
	s=$(echo "scale=$digits; $s0 + $ngamm" | bc)
	U $ngamm
	#echo "ngam=$ngamm tau0=$tau0, v0=$v0, s0=$s0, u=$u, v=$v, s=$s, y1max=$y1max, y1min=$y1min, delta=$delta"
	#echo "ystable=$ystable, yunstable=$yunstable, nystable=$nystable, nyunstable=$nyunstable, gammin=$gammini, gammfin=$gammfin, disty=$disty, ndisty=$ndisty"
	launchit
	auxystable=$(tail -n 1 stable_manifold.dat | awk '{print $2}')
	auxyunstable=$(tail -n 1 unstable_manifold.dat | awk '{print $2}')
	auxdisty=$(echo "scale=$digits; $auxystable- $auxyunstable" | bc)
	if [ $(echo "scale=$digits; $auxdisty*$disty<=0 "| bc) -eq 1 ]; then
		ystable=$nystable
		yunstable=$nyunstable
		disty=$ndisty
		gammini=$gammfin
	fi
	ndisty=$auxdisty
	andisty=$(abs $ndisty)
	echo "Error $andisty"
	gammfin=$ngamm
	nystable=$auxystable
	nyunstable=$auxyunstable
	echo $ngamm $tau0 $v0 $s0 0 $nystable $u $v $s $ndisty >> diff_stable_unstable.dat
	
#	newlimits_y1
	mv process_unstable process-unstable$i
	mv process_stable process-stable$i
done

echo "Finished"
