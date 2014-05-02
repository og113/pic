
reset
unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "phi(t,x)"
set xlabel "t"
set ylabel "x"
set zlabel "phi"
set grid
splot "./data/earlyphi_pic.dat" using 2:3:5 title 'phi(x,t)' with points
