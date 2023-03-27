set title 'Problem 13: part (a) and (b)'
set autoscale
set grid
set xlabel 'x'
m = 'myplot.dat'
set grid
plot m using 1:2 with lines title 'f(x)', m using 1:3 with linespoints title 'g(x)'