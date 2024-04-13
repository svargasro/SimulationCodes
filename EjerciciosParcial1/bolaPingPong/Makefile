all: bolaPingPong.x

bolaPingPong.x: bolaPingPong.cpp
	g++ bolaPingPong.cpp -o bolaPingPong.x
	./bolaPingPong.x | gnuplot


#	./bolaPingPong.x > data.txt
#	gnuplot plot.gp
clean:
	rm *.x *.txt
