IBMDIR = /home/kennedy/CPLEX_Studio221
CONCERTDIR = $(IBMDIR)/concert
CPLEXDIR = $(IBMDIR)/cplex
CPOPTDIR = $(IBMDIR)/cpoptimizer
GUROBIDIR = /home/kgaraujo/gurobi951/linux64

CFLAGS =	-O -DNDEBUG \
			-I$(CPOPTDIR)/include \
			-I$(CPLEXDIR)/include \
			-I$(CONCERTDIR)/include -fPIC -fstrict-aliasing -pedantic -std=c++14 -Wall -Wextra -fexceptions -frounding-math -Wno-long-long -m64

LDFLAGS =	-L$(CPOPTDIR)/lib/x86-64_linux/static_pic -lcp \
			-L$(CPLEXDIR)/lib/x86-64_linux/static_pic -lcplex -lilocplex \
			-L$(CONCERTDIR)/lib/x86-64_linux/static_pic \
			-lconcert  -lpthread -lm -ldl 

srcdir = /home/kgaraujo/Projetos/FJS

SOURCES =	Arc.cpp \
			GraphAlgorithms.cpp \
			Heuristics.cpp \
			Cenario.cpp \
			Modelos.cpp \
			SolutionGraph.cpp \
			SolutionNode.cpp \
			Tayebi.cpp \
			FJS.cpp \
			main.cpp
		
OBJECTS=$(SOURCES:.cpp=.o) 

TARGET = FJS


all: $(TARGET)


$(TARGET): main.cpp
	g++ -g -o $(TARGET) $(CFLAGS) $(SOURCES) $(LDFLAGS) 

clean:
	$(RM) $(TARGET) $(OBJECTS)