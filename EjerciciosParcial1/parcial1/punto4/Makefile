CXXFLAGS= -g -fsanitize=leak,address,undefined
all: parcial.x

parcial.x: parcial.cpp Makefile pPlot.py plot.gp
	g++ parcial.cpp -o parcial.x
	./parcial.x
	python3 pPlot.py
#	gnuplot plot.gp
#	open pPlot.png


#	./parcial.x > data.txt
#	gnuplot plot.gp


clean:
	rm *.x TablaydesvVSF_ConTorque.txt pPlot.png
