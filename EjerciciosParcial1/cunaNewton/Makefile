all: cunaNewton.x

cunaNewton.x: cunaNewton.cpp plot.gp initialPlot.py Makefile
	g++ cunaNewton.cpp -o cunaNewton.x
	./cunaNewton.x 2>"data.txt" 1| gnuplot
#	python3 initialPlot.py
	gnuplot plot.gp
	open initialPlot.png

clean:
	rm *.x *.txt
