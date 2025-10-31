#!/usr/bin/gnuplot

set terminal svg size 800,600
set output "plot.svg"
set xlabel "Magnetic Field Angle [deg]"
set ylabel "Particle Distance Change [μm]"
plot "TextData/TimeLog.csv" u (column("Phi")):(column("DeltaDistance")*1e6) t "Simulation"
