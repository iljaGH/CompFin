set logscale;
set xlabel 'level';
set ylabel 'relative error';
plot 'task9.dat' using 1:2 with lines title 'monte-carlo','task9.dat' using 1:3 with lines title 'trapezoidal','task9.dat' using 1:4 with lines title 'clenshaw-curtis','task9.dat' using 1:5 with lines title 'gauss-legendre'
