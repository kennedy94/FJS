# Directories for IBM CPLEX installation  
IBMDIR = /home/kennedy/CPLEX_Studio221
CONCERTDIR = $(IBMDIR)/concert
CPLEXDIR = $(IBMDIR)/cplex
CPOPTDIR = $(IBMDIR)/cpoptimizer

# Compilation flags  
CFLAGS =	-O -DNDEBUG \                 # Optimization and disable debugging  
			-I$(CPOPTDIR)/include \     # Include CP Optimizer header directories  
			-I$(CPLEXDIR)/include \     # Include CPLEX header directories  
			-I$(CONCERTDIR)/include \   # Include Concert Technology header directories  
			-fPIC -fstrict-aliasing -pedantic -std=c++14 -Wall -Wextra \  
			-fexceptions -frounding-math -Wno-long-long -m64  # Additional compilation options  

# Linking flags  
LDFLAGS =	-L$(CPOPTDIR)/lib/x86-64_linux/static_pic -lcp \  # Link CP Optimizer  
			-L$(CPLEXDIR)/lib/x86-64_linux/static_pic -lcplex -lilocplex \  # Link CPLEX  
			-L$(CONCERTDIR)/lib/x86-64_linux/static_pic \  # Link Concert Technology  
			-lconcert  -lpthread -lm -ldl  # Other required libraries  

# Directory where source files are located  
srcdir = /home/kgaraujo/Projetos/FJS

# List of source files in the project  
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
		
# List of object files to be generated from the sources  
OBJECTS=$(SOURCES:.cpp=.o) 

# Name of the final executable  
TARGET = FJS

# Main rule: compile the project  
all: $(TARGET)

# Rule to compile the executable  
$(TARGET): main.cpp
	g++ -g -o $(TARGET) $(CFLAGS) $(SOURCES) $(LDFLAGS) 

# Rule to clean compiled files  
clean:
	$(RM) $(TARGET) $(OBJECTS)
