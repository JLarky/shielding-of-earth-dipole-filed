#!/bin/bash

./3d.gnuplot
./2d.gnuplot
convert -density 110x110 -rotate 90 -crop 850x580+85+150 3d.ps 3d.png
