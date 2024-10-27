all:FyvsFm.png 




puntoF.x: puntoF.cpp
	g++ $^ -o $@
	rm -rf output/*
	parallel "./puntoF.x {} > ./output/w{#}.txt" ::: 100 200 300 400 500 600 700 800 900 1000

FyvsFm.png:  plot.py puntoF.x
	python3 plot.py

clean:
	rm -f *.o *.x *.dat
