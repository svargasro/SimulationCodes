CXXFLAGS= -g -fsanitize=leak,address,undefined
all: parcial.x

parcial.x: parcial.cpp Makefile pPlot.py pPlot.gp
	g++ parcial.cpp -o parcial.x
	./parcial.x 2>"data.txt" 1| gnuplot
	python3 pPlot.py
#gnuplot pPlot.gp
	open plot.png


#	./parcial.x > data.txt
#	gnuplot plot.gp


sanitizer:
	g++ ${CXXFLAGS} parcial.cpp -o parcial.x
	./parcial.x


clean:
	rm *.x *.txt
