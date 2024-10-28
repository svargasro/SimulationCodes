set terminal pdfcairo enhanced color font "Arial,10" size 5in, 4in
set output "Planet.pdf"

# Opciones de gráfico
set title "Planet orbit"
set xlabel "X-axis"
set ylabel "Y-axis"

# Datos o funciones
plot "data.txt" title "Planet orbit"

# Cerrar la salida
unset output
