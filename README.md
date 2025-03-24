Repository including the data and code used in the manuscripts https://arxiv.org/abs/2403.16766 and https://arxiv.org/abs/2403.16787. 

The code intends at solving the flexible job shop scheduling problem with sequencing flexibility and position-based learning effect.

Implementation of Contraint Programming em Mixed Integer Linear models, and heuristics/metaheuristic methods.

## Pre-requisites

IBM CPLEX 22.1

## Instance Format

Each instance consists of:
- The first line: number of operations, number of arcs, and number of machines.
- The next lines define arcs with their origin and destination.
- The remaining lines specify which machines can process each operation and their respective processing times.

### Example:

```plaintext
5 3 2  # 5 operations, 3 arcs, and 2 machines

0 1    # arc 0, origin and destination
0 2    # arc 1, origin and destination
3 4    # arc 2, origin and destination

2 0 1 1 2  # 2 machines can process operation 0, machine 0 with processing time 1, machine 1 with processing time 2
1 1 2      # 1 machine can process operation 1, machine 1 with processing time 2
2 0 1 1 2  # 2 machines can process operation 2, machine 0 with processing time 1, machine 1 with processing time 2
1 0 1      # 1 machine can process operation 3, machine 0 with processing time 1
2 0 1 1 2  # 2 machines can process operation 4, machine 0 with processing time 1, machine 1 with processing time 2
```
  
The precendence arcs can be shown like this:

  ![graph](https://github.com/user-attachments/assets/e97ff7c0-2016-4de1-bea9-175c26381aa8)

  
  
# How to run

Edit the Makefile to update the path to IBM CPLEX and the path to the source files.

Open terminal in the folder which the Makefile is

run the commands:

  chmod u+x makefile
  
  ./make

# Example
