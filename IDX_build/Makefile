CC=g++
CFLAGS=-O3 -m64  -std=c++0x  -fopenmp 
EXECUTABLE=MuGI_build
LIB=-ltbb
LDFLAGS=
SOURCES=data.cpp main.cpp bit_vectors.cpp  core.cpp  fasta.cpp  sa.cpp   sa2.cpp  sa_variants.cpp  suffix2.cpp  variant_list.cpp coding.cpp        sa1.cpp  sa3.cpp  suffix1.cpp      suffix3.cpp  vcf.cpp

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
