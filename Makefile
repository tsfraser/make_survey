# Uncomment to use GCC
#CXX = g++

# Uncomment to use Intel compiler
#CXX = icpc

# Uncomment to turn on aggressive optimizations
OPTFLAGS = -O3
GSL_INC= $(shell gsl-config --cflags) 
GSL_LINK= $(shell gsl-config --libs) 

CFLAGS= -Wall $(OPTFLAGS) -I./lib $(GSL_INC) 
CXXFLAGS = $(CFLAGS)
LDFLAGS = -lm $(GSL_LINK)

default: all 

all: make_survey

make_survey: make_survey.cpp lib/cuboid.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o
