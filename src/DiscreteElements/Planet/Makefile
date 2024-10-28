all: Planet.pdf

Planet.x: Planet.cpp
	g++ Planet.cpp -o Planet.x

data.txt: Planet.x
	./Planet.x > data.txt

Planet.pdf: data.txt plot.gp
	gnuplot plot.gp

clean:
	rm -f *.x *.pdf *.txt
