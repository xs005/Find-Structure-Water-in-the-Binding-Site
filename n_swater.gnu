reset
set terminal pngcairo
set grid
set xrange [*:]
set encoding iso_8859_1
set ylabel "Number of Structure Water"
set xlabel "Time (ns)"
set key bottom

set output "system_n_swater.png"
set title "Number of Structure Water in the Complex"
p '2.dat' u (column(1)/100):2 w l noti 
set output

