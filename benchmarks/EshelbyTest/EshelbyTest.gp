set terminal pngcairo enhanced size 1000, 800 font "Helvetica,15"
set datafile separator ","
set output 'result.png'

set style line 1 lt 3 lw 3 pt 17 lc rgb "red"
set style line 2 lt 3 lw 3 pt 6 lc rgb "green"
set style line 3 lt 3 lw 3 pt 18 lc rgb "blue"
set style line 4 lt 1 lw 3 pt 6 lc rgb "#fec287"

set xlabel "x/L_x"
set ylabel "Stress {/Symbol s} [MPa]"

plot "LineStresses.dat" u 1:2 ls 1 w linespoints title "{/Symbol s}_{r,OpenPhase}",\
     "LineStresses.dat" u 1:3 ls 2 w linespoints title "{/Symbol s}_{t,OpenPhase}",\
     "LineStresses.dat" u 1:4 ls 3 w linespoints title "{/Symbol s}_{r,Analytic}",\
     "LineStresses.dat" u 1:5 ls 4 w linespoints title "{/Symbol s}_{r,Analytic}",

quit
