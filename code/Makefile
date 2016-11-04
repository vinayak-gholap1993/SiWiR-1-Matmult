CXX = g++
CXXFLAGS = -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG

INCLUDES = 	
LDFLAGS =
LIBS =

# blas
INCLUDES += -I/usr/lib64/atlas/include/
LDFLAGS += -L/usr/lib64/atlas/
LIBS += -lcblas -latlas

# likwid
CXXFLAGS += -DUSE_LIKWID -pthread
INCLUDES += -I/usr/local/likwid/include/
LDFLAGS += -L/usr/local/likwid/lib/
LIBS += -llikwid

#TARGET = ex01_example
#OBJS = $(TARGET).o

HEADER = Timer.h Matrix.hpp IncludeFiles.hpp type.hpp
SRC =  Matrix
TEST = main

all: matmult

matmult: Makefile $(TEST).o $(SRC).o 
	$(CXX) $(CXXFLAGS) -o matmult $(TEST).o $(SRC).o  $(LDFLAGS) $(LIBS)

$(TEST).o: Makefile $(TEST).cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES)  $(TEST).cpp

$(SRC).o: Makefile $(SRC).cpp  $(HEADER)  Makefile 
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(SRC).cpp

clean:
	@$(RM) -rf *.o matmult
	@$(RM) -rf *.out
