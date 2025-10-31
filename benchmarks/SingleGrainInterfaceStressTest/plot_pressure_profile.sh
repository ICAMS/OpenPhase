#!/usr/bin/gnuplot

plot "tab_pressure_profile_al_al.ref" u "r_nm":"p_MPa" w lp title "Molecular Dynamics", "TextData/tab_data.csv" u (column("r")*1E+9):(column("p")*1E-6) w lp title "Phase-Field"
set title "Average radial pressure profile of an aluminium grain inside a aluminium matrix"
set xlabel "Radius   [nm]"
set ylabel "Pressure [Pa]"
set term svg
set output "fig_pressure_profile.svg"
replot
