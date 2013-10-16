# Makefile
#

#objects := $(patsubst %.cpp, %.o, $(wildcard *.cpp))
objects := HashSearch.o BlastStat.o Seg.o mergeUnit.o

CXX := g++

CXXFLAGs :=

ifeq ($(dbg), 1)
	CXXFLAGS += -g
else
	CXXFLAGS += -O3 -w
endif

ifeq ($(pg), 1)
	CXXFLAGS += -pg
else
	CXXFLAGS +=
endif

#INC := -I ~/include
INC := -I ./

#LIBS := -L ~/lib
LIBS := -L ./

all: rapsearch prerapsearch

rapsearch: main.cpp $(objects)
	$(CXX) $(CXXFLAGS) -o rapsearch main.cpp $(objects) $(INC) $(LIBS) -lboost_serialization -lpthread -lboost_thread -lboost_system -lboost_chrono -lrt

prerapsearch: BuildHash.cpp $(objects)
	$(CXX) $(CXXFLAGS) -o prerapsearch BuildHash.cpp $(objects) $(INC) $(LIBS) -lboost_serialization -lpthread -lboost_thread -lboost_system -lboost_chrono -lrt

$(objects) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)


clean:
	$(RM) *.o rapsearch prerapsearch
