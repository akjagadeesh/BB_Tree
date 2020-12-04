TARGET = benchmark
LIBS = -fopenmp -lm
CC = g++ 
CXXFLAGS = -mavx -fopenmp -std=c++11 -fno-tree-vectorize 
CFLAGS = -Wall -O3 -g -fopenmp -std=c++11 -mfpmath=both -ffast-math -fno-tree-vectorize -static-libstdc++

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
