all: Planeta.pdf


Planeta.x: Planeta.cpp
	g++ Planeta.cpp -o Planeta.x



# %.x: %.cpp
# 	g++ $^ -o $@
# 	Equivalente
#
#Si cambia sumupdown, entonces actualice data.txt


data.txt: Planeta.x
	./Planeta.x > data.txt


Planeta.pdf: data.txt plot.gp
	gnuplot plot.gp
