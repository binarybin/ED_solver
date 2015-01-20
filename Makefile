CC=g++
CFLAGS=-c -g -Wall -Wuninitialized -Wno-reorder -O3 -std=c++11
LDFLAGS=
LFLAGS=-lblas -llapack -llanczos -L/usr/local/gfortran/lib -lgfortran
SOURCES=hamiltonian.cpp hilbertspace.cpp interaction.cpp lanczos.cpp main.cpp measurement.cpp output.cpp solver.cpp diag_wrapper.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=ED_solver

all: $(SOURCES) $(EXECUTABLE) $(depend)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o
	rm $(EXECUTABLE)
