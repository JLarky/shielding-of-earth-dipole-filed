#!/usr/bin/gnuplot

set terminal png size 800,390
set xlabel 'GSM_X'
set ylabel 'GSM_Y'
unset zrange

set output '2d.png'
plot [20:-100][:] './points.dat' u 1:3 w d t 'X-Z plane'

set output 'XZ.png'
plot [20:-100][:] './XZ.dat' u 1:3 w p t 'X-Z plane, Y=0'
