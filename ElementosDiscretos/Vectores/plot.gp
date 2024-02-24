set term pdf

set out "Planeta.pdf"
set xlabel "x"
set ylabel "y"
plot "data.txt" with lines title "x vs. y"
