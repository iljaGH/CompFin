reset
set term png
set output "WS1T36.png"
set grid
set title '2d-plot Box-Muller method'
set xlabel 'z1'
set ylabel 'z2'
plot 'task6.dat' using 1:2
set xtic auto
set ytic auto