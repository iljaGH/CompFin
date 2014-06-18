set terminal postscript eps color;
set out '../LaTeX/task10Plot_10.eps';
set logscale;
set xlabel 'num nodes';
set ylabel 'relative error';
set grid;

plot 'task10_10.dat' using 1:2 with lines lt -1 lw 1 lc 0 title 'monte-carlo', \
'task10_10.dat' using 1:3 with lines lt -1 lw 1 lc 1 title 'trapezoidal', \
'task10_10.dat' using 1:4 with lines lt -1 lw 1 lc 2 title 'clenshaw-curtis', \
'task10_10.dat' using 1:5 with lines lt -1 lw 1 lc 3 title 'gauss-legendre'

set terminal postscript eps color;
set out '../LaTeX/task10Plot_0.eps';
set logscale;
set xlabel 'num nodes';
set ylabel 'relative error';

plot 'task10_0.dat' using 1:2 with lines lt -1 lw 1 lc 0 title 'monte-carlo', \
'task10_0.dat' using 1:3 with lines lt -1 lw 1 lc 1 title 'trapezoidal', \
'task10_0.dat' using 1:4 with lines lt -1 lw 1 lc 2 title 'clenshaw-curtis', \
'task10_0.dat' using 1:5 with lines lt -1 lw 1 lc 3 title 'gauss-legendre'
