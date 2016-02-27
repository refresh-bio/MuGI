CC=g++
CFLAGS=-O3 -m64  -std=c++0x -lrt   
EXECUTABLE=MuGI_search 
LIB=libaelf64.a
LDFLAGS=
SOURCES=SuffixArray.cpp  SuffixArray1.cpp  SuffixArray2.cpp  SuffixArray3.cpp  coding.cpp  main.cpp  seq_file.cpp  timer.cpp


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
