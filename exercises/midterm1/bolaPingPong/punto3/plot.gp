set term png
set out "initialPlot.png"
set xlabel 'x'
set ylabel 'y'


# Estilo de los puntos
set style line 1 pt 7 lc rgb "blue" # Puntos azules
set style line 2 pt 7 lc rgb "red"  # Puntos rojos

# Trazar los datos como puntos
plot 'data.txt' using 1:2 with points ls 1 title "pingpong", \
#     'data.txt' using 1:3 with points ls 2 title "raqueta"
