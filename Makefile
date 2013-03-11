# Uncomment to use GCC
#CXX = g++

# Uncomment to use Intel compiler
#CXX = icpc

# Uncomment to turn on aggressive optimizations
OPTFLAGS = -O3
GSL_INC= $(shell gsl-config --cflags)
GSL_LINK= $(shell gsl-config --libs)

CFLAGS= -Wall $(OPTFLAGS) -I./lib
CXXFLAGS = $(CFLAGS)
LDFLAGS = -lm

default: make_survey

all: make_survey tools

tools: mply_area mply_polyid mply_trim randbox

make_survey: make_survey.cpp lib/cuboid.cpp
	$(CXX) $(CXXFLAGS) $(GSL_INC) -o $@ $^ $(LDFLAGS) $(GSL_LINK)

mply_area: tools/mply_area.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

mply_polyid: tools/mply_polyid.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

mply_trim: tools/mply_trim.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

randbox: tools/randbox.c
	$(CC) $(CFLAGS) $(GSL_INC) -o $@ $^ $(LDFLAGS) $(GSL_LINK)

clean:
	rm -f *.o

bu-clean: clean
	rm -f *~ *.bak lib/*~ lib/*.bak tools/*~ tools/*.bak

real-clean: clean
	rm -f make_survey  mply_area mply_polyid mply_trim randbox
