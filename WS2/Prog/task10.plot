set term x11 0;
set logscale;
set xlabel 'level';
set ylabel 'relative error';
plot 'task10_0.dat' using 1:2 with lines title 'monte-carlo K=0','task10_0.dat' using 1:3 with lines title 'trapezoidal K=0','task10_0.dat' using 1:4 with lines title 'clenshaw-curtis K=0','task10_0.dat' using 1:5 with lines title 'gauss-legendre K=0'

set term x11 1;
set logscale;
set xlabel 'level';
set ylabel 'relative error';
plot 'task10_10.dat' using 1:2 with lines title 'monte-carlo  K=10','task10_10.dat' using 1:3 with lines title 'trapezoidal  K=10','task10_10.dat' using 1:4 with lines title 'clenshaw-curtis  K=10','task10_10.dat' using 1:5 with lines title 'gauss-legendre  K=10'
