all:cAvsRe.png Fyvst.png




puntoD.x: puntoD.cpp
	g++ $^ -o $@
	rm -rf output/*
	rm -rf output2/*
	parallel "./puntoD.x {}" ::: 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1

cAvsRe.png Fyvst.png:  plot.py plot2.py puntoD.x
	python3 plot.py
	python3 plot2.py

clean:
	rm -f *.o *.x *.dat
