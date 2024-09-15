GCC = g++
CFLAGS = --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = puntoE

all: $(TARGET).png

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(TARGET).cpp -o $(TARGET).out

FxFy.dat: $(TARGET).out
	rm -f *.dat
	./$(TARGET).out 32 2>> "FxFy.dat"
	./$(TARGET).out 28 2>> "FxFy.dat"
	./$(TARGET).out 24 2>> "FxFy.dat"
	./$(TARGET).out 16 2>> "FxFy.dat"
	./$(TARGET).out 12 2>> "FxFy.dat"


$(TARGET).png: FxFy.dat plotFyFx.py
	python3 plotFyFx.py

.PHONY: clean

clean:
	rm *.out *.dat *.png
