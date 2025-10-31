Magnetoactive Elastomer Linear
============================
This is a benchmark for the coupled magneto-mechanical boundary value Problem
in magneto-active elastomers. Two magnetizable particles are embedded in a
polymer matrix. The particles are magnetized by an external magnetic field. The
resulting displacement of the particles is calculated as function function of
the direction of the external magnetic field.

The benchmark results are written into TextData/TimeLog.csv and be compared to
the results published [1].

The results can be plotted using gnuplot on Linux with
./plot.p
which produced a svg plot plot.svg.

How to Execute:
---------------
./MagnetoactiveElastomerLinear

TODOs:
---------------
+ Add compare.sh script for automated comparing of results

Authors:
-------
raphael.schiedung@rub.de

References:
-------
[1] Metsch, Philipp, et al. "Benchmark for the coupled magneto-mechanical
boundary value problem in magneto-active elastomers." Materials 14.9 (2021): 2380.
