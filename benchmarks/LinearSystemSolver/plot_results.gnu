#!/usr/bin/gnuplot -persist
set datafile separator ','
set xlabel 'x'
set ylabel 'y'
plot "./TextData/Gauss.csv"             u 1:(column('Gauss')),\
     "./TextData/Jacobi.csv"            u 1:(column('Jacobi')),\
     "./TextData/GaussSeidel.csv"       u 1:(column('GaussSeidel')),\
     "./TextData/GradientDecent.csv"    u 1:(column('GradientDecent')),\
     "./TextData/ConjugateGradient.csv" u 1:(column( 'ConjugateGradient')),\
     "./TextData/BiconjugateGradient.csv" u 1:(column('BiconjugateGradient')),\
                                       "" u 1:(column('Analytic')) w l


