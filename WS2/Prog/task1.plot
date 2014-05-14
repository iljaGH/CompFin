reset
set terminal epslatex color
set out 'task1Plot.tex'
set xlabel "$\\sigma$";
set ylabel "mean-estimate";
set grid
set xrange [-0.1:0.9]
unset key
plot 'task1.dat' using 1:2 pt 7 ps 2
set out
