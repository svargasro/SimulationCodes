GCC = g++
CFLAGS = --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = puntoAD

all: $(TARGET).png

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(TARGET).cpp -o $(TARGET).out

wind.dat: $(TARGET).out
	./$(TARGET).out

$(TARGET).png FX_FY_t.png: wind.dat plot.py plotFyFx.py
	python3 plot.py
	python3 plotFyFx.py

.PHONY: clean

clean:
	rm *.out *.dat
