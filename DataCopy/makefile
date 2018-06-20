CXXFLAGS=-std=c++11

all:NDCopy NDCopy2 main.cpp
	g++ main.cpp $(CXXFLAGS) NDCopy.o NDCopy2.o -o exe


NDCopy:NDCopy.cpp
	g++ NDCopy.cpp $(CXXFLAGS) -c

NDCopy2:NDCopy2.cpp
	g++ NDCopy2.cpp $(CXXFLAGS) -c


clean:
	rm *.o
