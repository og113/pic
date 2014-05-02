
reset
unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "action(a): epsilon = 0.0001, R = 100, box = 200"
set xlabel "a"
set ylabel "action"
set grid
set key below
f(x) = a*x + b
a=-1.0; b=1.0;
fit f(x) "./data/action_pic1.dat" using (100/$2):($7) via a, b
plot "./data/action_pic1.dat" using (100/$2):($7) title 'lattice result' with points, \
	f(x) title sprintf('fit: f(x) = %.2f*x + %.2f',a,b)
