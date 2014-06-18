set terminal postscript eps color;
set out '../LaTeX/task4Plot.eps';
set xlabel 'N';
set ylabel "approximation to fair value";
set logscale x;
set logscale y;
set grid
plot 'task4_1.dat' using 1:2 with lines title 'run 1' lt -1 lw 1 lc 0, 'task4_2.dat' using 1:2 with lines title 'run 2' lt -1 lw 1 lc 1, 'task4_3.dat' using 1:2 with lines title 'run 3' lt -1 lw 1 lc 2, 'task4_4.dat' using 1:2 with lines title 'run 4' lt -1 lw 1 lc 3, 'task4_5.dat' using 1:2 with lines title 'run 5' lt -1 lw 1 lc 4;

set terminal postscript eps color;
set out '../LaTeX/task4Plot_err.eps';
set xlabel 'N';
set ylabel "relative error";
set logscale x;
set logscale y;
set grid
plot 'task4_1_conv.dat' using 1:2 with lines title 'run 1' lt -1 lw 1 lc 0, 'task4_2_conv.dat' using 1:2 with lines title 'run 2' lt -1 lw 1 lc 1, 'task4_3_conv.dat' using 1:2 with lines title 'run 3' lt -1 lw 1 lc 2, 'task4_4_conv.dat' using 1:2 with lines title 'run 4' lt -1 lw 1 lc 3, 'task4_5_conv.dat' using 1:2 with lines title 'run 5' lt -1 lw 1 lc 4;
