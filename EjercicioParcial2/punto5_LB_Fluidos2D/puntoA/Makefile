all: plotA.pdf

%.x: %.o
	g++ $^ -o $@

wind.dat: puntoA.x
	./puntoA.x

plotA.pdf: wind.dat plot.py
	python3 plot.py

clean:
	rm -f *.o *.x *.dat
