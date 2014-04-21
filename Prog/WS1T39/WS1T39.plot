#set term x11 0;
set logscale xy;
set ylabel "estimation error";
set title "estimation of mu";
plot 'task9_1.dat' using 1:1 title columnheader(1) with lines linecolor rgb "white",'task9_1.dat' using 1:2 with lines title columnheader(2),'task9_2.dat' using 1:2 with lines title columnheader(2), 'task9_3.dat' using 1:2 with lines title columnheader(2), x**(-0.5) title "N^(-0.5)";

#set term x11 1;
set logscale xy;
set ylabel "estimation error";
set title "estimation of sigma given mu";
plot 'task9_1.dat' using 1:1 title columnheader(1) with lines linecolor rgb "white",'task9_1.dat' using 1:3 with lines title columnheader(2),'task9_2.dat' using 1:3 with lines title columnheader(2), 'task9_3.dat' using 1:3 with lines title columnheader(2), x**(-0.5) title "N^(-0.5)";
