CXX = g++
CXXFLAGS = -O3 -m64 -std=c++17 -fopenmp

SOURCE = sdis-*.cpp

all: rss

bin:
	mkdir -p bin

rss: rss-count

rss-count: bin
	$(CXX) $(CXXFLAGS) -o bin/$@ $@.cpp $(SOURCE)

clean:
	rm -rf bin
