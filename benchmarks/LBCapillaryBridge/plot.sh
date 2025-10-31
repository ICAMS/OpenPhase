#!/usr/bin/gnuplot
#set terminal postscript eps enhanced color font 'helvetica,12'
set terminal pdf
set output 'plot.pdf'
set datafile separator ','

set title 'Capillarity bridge - Wall distance'
set xlabel 'Rescaled Time'
set ylabel 'Rescaled Wall distance'
set yrange [1:2]
set xrange [0:0.8]
plot "TextData/Log_Observables.csv" u (column("TimeRescaled")):(column("WallGapRescaled")) title 'Simulation', \
    ""u (column("TimeRescaled")):(column("WallGapRescaledEx")) title 'Theory' w l, \
