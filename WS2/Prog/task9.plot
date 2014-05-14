reset
set terminal epslatex color
set out 'task9Plot.tex'
set logscale;
set xlabel 'num. nodes';
set ylabel 'relative error';

set style line 1 lc rgb '#8b1a0e' pt 6 ps 1 lt 1 lw 2 
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2
set style line 3 lc rgb 'blue' pt 6 ps 1 lt 1 lw 2
set style line 4 lc rgb 'black' pt 6 ps 1 lt 1 lw 2

set style line 11 lc rgb '#ffffff' lt 1
set border 0 back ls 11
#set tics out nomirror scale 0,0.001
#set format ''

set mxtics
set mytics
set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 12, ls 13

plot 'task9.dat' using 1:2 with lp ls 1 title 'monte-carlo', \
'task9.dat' using 1:3 with lp ls 2 title 'trapezoidal', \
'task9.dat' using 1:4 with lp ls 3 title 'clenshaw-curtis', \
'task9.dat' using 1:5 with lp ls 4 title 'gauss-legendre'
set out