Repository including the data and code used in the manuscripts https://arxiv.org/abs/2403.16766 and https://arxiv.org/abs/2403.16787. 

The code intends at solving the flexible job shop scheduling problem with sequencing flexibility and position-based learning effect.

Implementation of Contraint Programming em Mixed Integer Linear models, and heuristics/metaheuristic methods.

# Pre-requisites

IBM CPLEX 22.1

# Instance format

  \#operations \#arcs \#machines
  
  for each arc:
    origin destination
    
  for each operation:  
    \#machines that can process operation
    for each machine that can process operation:
      machine processing-time 
  
  Example:
  
  5 3 2
  
  0 1
  
  0 2
  
  3 4

  2 0 1 1 2
  
  1 1 2
  
  2 0 1 1 2
  
  1 0 1
  
  2 0 1 1 2
  
  
![graph](https://github.com/user-attachments/assets/e97ff7c0-2016-4de1-bea9-175c26381aa8)

  
  
# How to run

Edit the Makefile to update the path to IBM CPLEX and the path to the source files.

Open terminal in the folder which the Makefile is

run the commands:

  chmod u+x makefile
  
  ./make

# Example
