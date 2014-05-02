
reset
unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "x,v"
set xlabel "t"
set ylabel "x,v"
set grid
plot "./data/lf.dat" using 1:2 title 'x(t)' with points, \
	"./data/lf.dat" using 1:3 title 'v(t)' with points
