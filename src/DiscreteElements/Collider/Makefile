SRC = Collider.cpp
EXE = $(SRC:.cpp=.x)

all: $(EXE)

# Rule to compile the executable
$(EXE): $(SRC)
	g++ $(SRC) -o $(EXE)
	./$(EXE) | gnuplot

clean:
	rm -f *.x
