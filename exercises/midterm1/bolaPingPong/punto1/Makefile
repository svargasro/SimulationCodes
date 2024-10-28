CXXFLAGS= -g -fsanitize=leak,address,undefined
all: bolaPingPong.x

bolaPingPong.x: bolaPingPong.cpp
	g++ bolaPingPong.cpp -o bolaPingPong.x
	./bolaPingPong.x 2>"data.txt" 1| gnuplot
	python3 initialPlot.py
#gnuplot plot.gp
	open initialPlot.png



#	./bolaPingPong.x > data.txt
#	gnuplot plot.gp


sanitizer:
	g++ ${CXXFLAGS} bolaPingPong.cpp -o bolaPingPong.x
	./bolaPingPong.x


clean:
	rm *.x *.txt
