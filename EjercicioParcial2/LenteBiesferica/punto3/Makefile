GCC = g++
CFLAGS = --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = lenteBiesferica

all: $(TARGET).png

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(TARGET).cpp -o $(TARGET).out

Waves2D.dat: $(TARGET).out
	./$(TARGET).out

$(TARGET).png: Waves2D.dat
	gnuplot plotWaves.gp

.PHONY: clean

clean:
	rm *.out *.dat *.png
