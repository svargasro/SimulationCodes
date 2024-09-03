# Compilador
CXX = g++



# Obtiene una lista de todos los archivos .cpp
SRC = $(wildcard *.cpp)

# Convierte los nombres de archivos .cpp a nombres de ejecutables .x
EXEC = $(SRC:.cpp=.x)

# Regla por defecto que compila todos los archivos .cpp
all: $(EXEC)

# Regla para convertir .cpp a .x
%.x: %.cpp
	$(CXX) -o $@ $<
	./$@

# Limpiar archivos generados
clean:
	rm -f $(EXEC)
