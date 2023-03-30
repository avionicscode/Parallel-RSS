CXX = g++
CXXFLAGS = -O3 -m64 -fopenmp -std=c++17 

SOURCE = sdis-*.cpp

all: rss

bin:
	mkdir -p bin

rss: rss-time

rss-time: bin
	$(CXX) $(CXXFLAGS) -o bin/$@ $@.cpp $(SOURCE)

clean:
	rm -rf bin
