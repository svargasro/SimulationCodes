SRC = LB_Fluids_D2Q9.cpp
EXE = $(SRC:.cpp=.x)
PLOT = plotFluids.gp
DATA = WindChannel.dat

all: $(DATA)

# Rule to compile the executable
$(EXE): $(SRC)
	g++ $(SRC) -o $(EXE)
	./$(EXE)

$(DATA): $(EXE) $(PLOT)
	gnuplot $(PLOT)
clean:
	rm -f *.x *.png *.dat
