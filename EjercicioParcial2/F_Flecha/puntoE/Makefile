GCC = g++
CFLAGS = --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = puntoE

all: $(TARGET).png

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(TARGET).cpp -o $(TARGET).out

Fxb.dat: $(TARGET).out
	rm -f *.dat
	./$(TARGET).out 32 2>> "Fxb.dat"
	./$(TARGET).out 28 2>> "Fxb.dat"
	./$(TARGET).out 24 2>> "Fxb.dat"
	./$(TARGET).out 16 2>> "Fxb.dat"
	./$(TARGET).out 12 2>> "Fxb.dat"



$(TARGET).png: Fxb.dat plotFxb.py
	python3 plotFxb.py

.PHONY: clean

clean:
	rm *.out *.dat
