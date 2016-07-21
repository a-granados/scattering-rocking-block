#!/bin/bash

###This is just to try this function before using it in find_het4.sh

function newlimits2_y1(){
###This is still not working.


##In this version we find in diff_stable_unstable.dat the interval
##where ngam is located. Then we look for y1min and y1max in unstable_manifold.dat
##and stable_manifold.dat
##We should have calculated at least two iterates.

local numdata
local file1

local k
local aux


file1=ordered_gammas2

numdata=$(cat $file1.dat | wc -l)

aux=$(tail -n +1 $file1.dat | head -n 1 | awk '{print $1}')

if [ $(echo "scale=$digits; $ngamm<$aux" | bc) -eq 1 ]; then
	echo "$ngamm $i" > $file1.tmp
	y1max=$firsty1max
	y1min=$firsty1min
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
		y1max=$firsty1max
		y1min=$firsty1min
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
		aux2=$(max $y1 $y2)
		aux3=$(max $y3 $y4)
		y1max=$(max $aux2 $aux3)
		aux2=$(min $y1 $y2)
		aux3=$(min $y3 $y4)
		y1min=$(min $aux2 $aux3)
		head -n $k $file1.dat > $file1.tmp2
		let aux=$numdata-$k
		tail -n $aux $file1.dat > $file1.tmp3
		cat $file1.tmp2 $file1.tmp $file1.tmp3 > $file1.dat
		rm -f $file1.tmp*
	fi
fi
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

firsty1max=1.1
firsty1min=0.9

digits=60

i=3
ngamm=.825287809941879781139948820394827823697781341575894138433948

newlimits2_y1

echo $y1min
echo $y1max




