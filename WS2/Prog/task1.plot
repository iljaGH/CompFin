reset;
set terminal postscript eps color;
set out '../LaTeX/task1Plot.eps';
set xlabel "sigma";
set ylabel "mean estimate";
set grid
set xrange [-0.1:0.9]
plot 'task1.dat' using 1:2 pt 7 ps 2 notitle
