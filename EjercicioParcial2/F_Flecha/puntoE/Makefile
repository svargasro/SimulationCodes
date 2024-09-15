GCC = g++
CFLAGS = --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = puntoE

VALUES = $(shell seq 12 4 40)

all: $(TARGET).png

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(TARGET).cpp -o $(TARGET).out

Fxb.dat: $(TARGET).out
	rm -f *.dat
	@for value in $(VALUES); do \
	echo "Ejecutando con b=$$value"; \
	./$(TARGET).out $$value 2>> "Fxb.dat"; \
	done



$(TARGET).png: Fxb.dat plotFxb.py
	python3 plotFxb.py

.PHONY: clean

clean:
	rm *.out *.dat
