SRC = Waves_D2Q5.cpp
EXE = $(SRC:.cpp=.x)
PLOT = plotWaves.gp
DATA = Waves2D.dat

all: $(DATA)

# Rule to compile the executable
$(EXE): $(SRC)
	g++ $(SRC) -o $(EXE)
	./$(EXE)

$(DATA): $(EXE)
	gnuplot $(PLOT)
clean:
	rm -f *.x *.png *.dat
