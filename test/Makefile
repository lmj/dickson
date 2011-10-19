
CXX = g++

OPTIMIZE = -O3 -ffast-math -fno-exceptions -fno-rtti -fomit-frame-pointer -DNDEBUG -s

INC_DIRS = -I../dickson

CFLAGS = $(INC_DIRS) $(OPTIMIZE) $(WARNINGS) $(EXTRA_CFLAGS)

all : laws basic

laws : laws.cxx common.hxx ../dickson/dickson.hxx
	$(CXX) $(CFLAGS) laws.cxx -o laws

basic : basic.cxx common.hxx ../dickson/dickson.hxx
	$(CXX) $(CFLAGS) basic.cxx -o basic

clean : 
	rm -f laws basic

.PHONY : all clean

