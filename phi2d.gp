
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
set key below
splot '<paste ./data/phi_pic0.dat ./data/phi0.dat' using 2:3:($5-$9) title 'difference' with points lt rgb "blue"
