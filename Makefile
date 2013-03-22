# Uncomment to use GCC
#CXX = g++

# Uncomment to use Intel compiler
#CXX = icpc

OPTFLAGS = -O3

# edit these if gsl-config is not in your default path
GSL_INC= $(shell gsl-config --cflags)
GSL_LINK= $(shell gsl-config --libs)

CFLAGS= -Wall $(OPTFLAGS) -I./lib
CXXFLAGS = $(CFLAGS)
LDFLAGS = -lm

default: make_survey

all: make_survey tools

tools: check_survey mply_area mply_polyid mply_trim randbox trim_by_index

make_survey: src/make_survey.cpp lib/cuboid.cpp
	$(CXX) $(CXXFLAGS) $(GSL_INC) -o $@ $^ $(LDFLAGS) $(GSL_LINK)

check_survey: src/check_survey.cpp lib/cuboid.cpp
	$(CXX) $(CXXFLAGS) $(GSL_INC) -o $@ $^ $(LDFLAGS) $(GSL_LINK)

mply_area: src/mply_area.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

mply_polyid: src/mply_polyid.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

mply_trim: src/mply_trim.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

randbox: src/randbox.c
	$(CC) $(CFLAGS) $(GSL_INC) -o $@ $^ $(LDFLAGS) $(GSL_LINK)

trim_by_index: src/trim_by_index.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o

bu-clean: clean
	rm -f src/*~ src/*.bak lib/*~ lib/*.bak 

real-clean: clean
	rm -f make_survey  check_survey mply_area mply_polyid mply_trim randbox trim_by_index
