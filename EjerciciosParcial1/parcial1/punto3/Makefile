CXXFLAGS= -g -fsanitize=leak,address,undefined
all: parcial.x

parcial.x: parcial.cpp Makefile pPlot.py plot.gp
	g++ parcial.cpp -o parcial.x
	./parcial.x
	python3 pPlot.py
	python3 plotFrecuencia2.py
#	gnuplot plot.gp
#	open pPlot.png


#	./parcial.x > data.txt
#	gnuplot plot.gp


sanitizer:
	g++ ${CXXFLAGS} parcial.cpp -o parcial.x
	./parcial.x


clean:
	rm *.x *.txt
