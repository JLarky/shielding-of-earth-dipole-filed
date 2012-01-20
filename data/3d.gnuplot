#!/usr/bin/gnuplot

set term postscript enhanced color #portrait
set output "3d.ps"
set view 80,105,1,1
set size .9,1.1

set zrange [-29:30]
set xlabel 'GSM_X' #offset -1,9
set ylabel 'GSM_Y' offset 2,-1
set zlabel 'GSM_Z' offset -1,9
splot [:][-33:33] './points.dat' w d t ''
