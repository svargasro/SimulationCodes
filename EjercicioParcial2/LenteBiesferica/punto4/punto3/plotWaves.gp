set terminal jpeg
set output "lenteBiesferica.jpg"
set pm3d map
set size ratio 1
set xlabel "X" font "Arial,14"
set ylabel "Y" font "Arial,14"
set autoscale xfix
set autoscale yfix
splot "Waves2D.dat"


# set term png
# set out "WavesPlot.png"
# set xlabel "X"
# set ylabel "Y"
# set zlabel "Z"
# splot "Waves2D.dat" w l
