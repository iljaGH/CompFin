set logscale;
set xlabel 'N';
set ylabel 'm_n';
plot 'task4_1.dat' using 1:2 with lines title 'run 1', 'task4_2.dat' using 1:2 with lines title 'run 2', 'task4_3.dat' using 1:2 with lines title 'run 3', 'task4_4.dat' using 1:2 with lines title 'run 4', 'task4_5.dat' using 1:2 with lines title 'run 5';
