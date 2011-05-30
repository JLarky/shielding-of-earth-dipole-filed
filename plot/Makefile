all: run

PROGNAME=plot
include ../inc/makefile.in


lines: lines.ps

lines.ps: lines.dat plot.gnuplot
	./plot.gnuplot
