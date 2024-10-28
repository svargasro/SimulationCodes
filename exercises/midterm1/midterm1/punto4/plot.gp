set term png
set out "pPlot.png"
set xlabel 'x'
set ylabel 'y'


plot 'data.txt' w p
