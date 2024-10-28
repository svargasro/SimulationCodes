set terminal jpeg
set output "lenteBiesferica.jpg"
set pm3d map
set size ratio -1
set xlabel "X" font "Arial,14"
set ylabel "Y" font "Arial,14"
set autoscale xfix
set autoscale yfix
#set xrange [32:96]  # Define el rango del eje X
#set xtics 1         # Establece los tics en intervalos de 1
set xtics format "" # Oculta los valores numéricos en el eje X


splot "Waves2D.dat"


# set term png
# set out "WavesPlot.png"
# set xlabel "X"
# set ylabel "Y"
# set zlabel "Z"
# splot "Waves2D.dat" w l
