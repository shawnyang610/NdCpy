CXX=g++-8
CXXFLAGS=-std=c++11 -g

all:NDCopy NDCopy2 main.cpp
	$(CXX) main.cpp $(CXXFLAGS) -c
	$(CXX) $(CXXFLAGS) *.o -o exe


NDCopy:NDCopy.cpp
	$(CXX) NDCopy.cpp $(CXXFLAGS) -c

NDCopy2:NDCopy2.cpp
	$(CXX) NDCopy2.cpp $(CXXFLAGS) -c


clean:
	rm *.o
