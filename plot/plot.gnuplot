#!/usr/bin/gnuplot

set term postscript enhanced color portrait
set output "lines.ps"
set size 2.5,1.3

#plot [10:-10][-10:10] './lines.dat' u 1:3 w p

set size 1.1,0.5
plot [20:-90][-40:40] './lines.dat' u 1:3 w d t 'field lines', '../data/XZ.dat' u 1:3 t 'magnetopause'
#plot [30:-30][-20:20] './lines.dat' u 1:3 w d t 'fiels lines', '../data/XZ.dat' u 1:3
