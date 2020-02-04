# This make file compiles the C/C++ executables used by examples/Example5.java
# and examples/Example6.java.
#
# This make file assumes it is being run on a Unix-like environment with GCC
# and Make installed.  Windows users can install MinGW/MSYS or Cygwin.
#
# Run 'make' to build all files

CXX = mpiCC
CXXFLAGS = -c -Wall -g
INC = -I /opt/boost/1.54.0/include
LIBS = -lm
UNAME_S = $(shell uname -s)

ifneq (, $(findstring SunOS, $(UNAME_S)))
    LIBS += -lnsl -lsocket -lresolv
endif

all: portfolio

portfolio: portfolio.o moeaframework.o 
	$(CXX) portfolio.o moeaframework.o -o portfolio

moeaframework.o: ./moeaframework.c ./moeaframework.h
	$(CXX) $(CXXFLAGS) ./moeaframework.c

portfolio.o: portfolio.cpp ./moeaframework.h ./boostutil.h 
	$(CXX) $(CXXFLAGS) $(INC) portfolio.cpp -o portfolio.o
 


