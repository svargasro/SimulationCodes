all: Planeta.pdf


PlanetaConstructor.x: PlanetaConstructor.cpp
	g++ PlanetaConstructor.cpp -o PlanetaConstructor.x



# %.x: %.cpp
# 	g++ $^ -o $@
# 	Equivalente
#
#Si cambia sumupdown, entonces actualice data.txt


data.txt: PlanetaConstructor.x
	./PlanetaConstructor.x > data.txt


Planeta.pdf: data.txt plot.gp
	gnuplot plot.gp
