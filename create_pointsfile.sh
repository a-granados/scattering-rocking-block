#!/bin/bash

digits=60

x1=$(echo $1 |  sed 's/^\-\./-0\./')
x2=$(echo $2 |  sed 's/^\-\./-0\./')
y2=$(echo $3 |  sed 's/^-\./-0\./')
t0=$(echo $4 |  sed 's/^\-\./-0\./')
y1min=$(echo $5 |  sed 's/^\-./-0\./')
y1max=$(echo $6 |  sed 's/^\-\./-0\./')
ypoints=$7
nfile=$8

if [ $ypoints -le 1 ]; then
	echo "Increase number of points in create_pointsfile.sh"
	exit
fi
hy1=$(echo "scale=$digits; ($y1max - $y1min)/($ypoints-1)" | bc)

y1=$y1min
for i in `seq 1 1 $ypoints`; do
#	echo $x1 $y1 $x2 $y2 $t0 >> points.dat
	y1=$(echo "scale=$digits; $y1min+($i-1)*$hy1" | bc | sed 's/^\-\./-0\./'  )
	echo $x1 $y1 $x2 $y2 $t0 >> $nfile
done

