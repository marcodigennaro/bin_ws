#!/bin/sh

sumX=0
sumX2=0
sumY=0
sumXY=0
count=0
while IFS=" " read -r x y 
do
    sumX=`awk "BEGIN {print $sumX+$x}"`
    sumY=`awk "BEGIN {print $sumY+$y}"`
    sumX2=`awk "BEGIN {print $sumX2+$x*$x}"`
    sumXY=`awk "BEGIN {print $sumXY+$x*$y}"`
    count=$(expr $count + 1)
done <  $1

b=`awk "BEGIN {print ($count*$sumXY-$sumX*$sumY)/($count*$sumX2 - $sumX*$sumX) }"`
a=`awk "BEGIN {print ($sumY - $b*$sumX)/$count }"`
mean=`awk "BEGIN {print ($sumY/$count)}"`
d_square=0

while IFS=" " read -r x y 
do
    d_square=`awk "BEGIN {print ($y - $mean)*($y - $mean)}"`
done < $1
std=`awk "BEGIN { print sqrt($d_square/$count)}"`

echo "a=$a, b=$b, mean=$mean, std=$std"

