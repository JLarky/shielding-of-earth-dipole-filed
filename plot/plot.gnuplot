#!/usr/bin/gnuplot

set term postscript enhanced color portrait
set output "lines.ps"
set size 2.5,1.3

#plot [10:-10][-10:10] './lines.dat' u 1:3 w p

set size 1.1,0.5
plot [40:-90][-60:60] './lines.dat' u 1:3 w d t 'dipole+cusp', '../data/XZ.dat' u 1:3
plot [30:-30][-20:20] './lines.dat' u 1:3 w d t 'dipole+cusp', '../data/XZ.dat' u 1:3
