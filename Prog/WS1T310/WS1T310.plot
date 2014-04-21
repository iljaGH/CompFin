#set term x11 0;
plot 'task10_t1.dat' using 1:2 with lines title 'delta t=0.5','task10_t1.dat' using 1:3 with lines title 'delta t=0.5','task10_t1.dat' using 1:4 with lines title 'delta t=0.5','task10_t2.dat' using 1:2 with lines title 'delta t=0.01','task10_t2.dat' using 1:3 with lines title 'delta t=0.01','task10_t2.dat' using 1:4 with lines title 'delta t=0.01';

#set term x11 1;
plot x;
plot 'task10_t1.dat' using 1:5 with lines title 'delta t=0.5','task10_t1.dat' using 1:6 with lines title 'delta t=0.5','task10_t1.dat' using 1:7 with lines title 'delta t=0.5','task10_t2.dat' using 1:5 with lines title 'delta t=0.01','task10_t2.dat' using 1:6 with lines title 'delta t=0.01','task10_t2.dat' using 1:7 with lines title 'delta t=0.01';