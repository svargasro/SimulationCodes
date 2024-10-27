all: plotE.pdf

%.x: %.o
	g++ $^ -o $@

wind.dat: puntoE.x
	./puntoE.x

plotE.pdf: wind.dat plot.py
	python3 plot.py

clean:
	rm -f *.o *.x *.dat