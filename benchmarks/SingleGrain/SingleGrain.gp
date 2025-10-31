set terminal pngcairo size 800, 800 font "Helvetica,14"
set datafile separator ","
set output 'R2.png'

set style line 1 lt 1 lw 3 pt 3 lc rgb "red"
set style line 2 lt 1 lw 3 pt 3 lc rgb "blue"
set style line 3 lt 1 lw 3 pt 3 lc rgb "green"
set style line 4 lt 1 lw 3 pt 3 lc rgb "black"

set xlabel "R^2"
set ylabel "Timesteps"

plot "R_2_graph.dat" u 1:2 ls 1 w l title "Analytic", \
     "R_2_graph.dat" u 1:3 ls 2 w l title "OpenPhase"
quit
