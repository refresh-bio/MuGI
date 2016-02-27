CC=g++
CFLAGS=-O3 -m64  -std=c++0x  
EXECUTABLE=MuGI_genEx
LIB=
LDFLAGS=
SOURCES=main.cpp coding.cpp  
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) $(LIB) -o $@
	cp $(EXECUTABLE) ../

.cpp.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm $(OBJECTS)
	rm $(EXECUTABLE)
	rm ../$(EXECUTABLE)
