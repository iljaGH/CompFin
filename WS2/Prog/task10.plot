set terminal postscript eps color;
set out '../LaTeX/task10Plot_10.eps';
set logscale;
set xlabel 'N nodes';
set ylabel 'relative error';
set grid;

plot 'task10_10.dat' using 1:2 with lines lt -1 lw 1 lc 0 title 'Monte Carlo', \
'task10_10.dat' using 1:3 with lines lt -1 lw 1 lc 1 title 'Trapezoidal', \
'task10_10.dat' using 1:4 with lines lt -1 lw 1 lc 2 title 'Clenshaw-Curtis', \
'task10_10.dat' using 1:5 with lines lt -1 lw 1 lc 3 title 'Gauss-Legendre'

set terminal postscript eps color;
set out '../LaTeX/task10Plot_0.eps';
set logscale;
set xlabel 'N nodes';
set ylabel 'relative error';

plot 'task10_0.dat' using 1:2 with lines lt -1 lw 1 lc 0 title 'Monte Carlo', \
'task10_0.dat' using 1:3 with lines lt -1 lw 1 lc 1 title 'Trapezoidal', \
'task10_0.dat' using 1:4 with lines lt -1 lw 1 lc 2 title 'Clenshaw-Curtis', \
'task10_0.dat' using 1:5 with lines lt -1 lw 1 lc 3 title 'Gauss-Legendre'
