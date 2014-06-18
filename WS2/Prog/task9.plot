reset
set terminal postscript eps color;
set out '../LaTeX/task9Plot.eps';
set logscale;
set xlabel 'num. nodes';
set ylabel 'relative error';

set grid
plot 'task9.dat' using 1:2 with lines title 'monte-carlo' lt -1 lw 1 lc 0, \
'task9.dat' using 1:3 with lines title 'trapezoidal' lt -1 lw 1 lc 1, \
'task9.dat' using 1:4 with lines title 'clenshaw-curtis' lt -1 lw 1 lc 2, \
'task9.dat' using 1:5 with lines title 'gauss-legendre' lt -1 lw 1 lc 3
set out
