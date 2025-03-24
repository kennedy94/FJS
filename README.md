# Flexible Job Shop Scheduling Problem  

This repository contains the data and code used in the manuscripts:  
- [arXiv:2403.16766](https://arxiv.org/abs/2403.16766)  
- [arXiv:2403.16787](https://arxiv.org/abs/2403.16787)  

The code aims to solve the **Flexible Job Shop Scheduling Problem** with sequencing flexibility and position-based learning effects.  
It includes implementations of **Constraint Programming**, **Mixed Integer Linear Programming (MILP)** models, and **heuristic/metaheuristic** methods.  

---

## 📌 Prerequisites  

- **IBM CPLEX 22.1**  

---

## 📄 Instance Format  

Each instance consists of:  
- The **first line**: Number of operations, number of arcs, and number of machines.  
- The **next lines**: Arcs defined by their origin and destination.  
- The **remaining lines**: Machines that can process each operation and their respective processing times.  

### 🔹 Example:  

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
### 🔹 Precedence arcs representation

The precedence arcs can be visualized as follows:

  ![graph](https://github.com/user-attachments/assets/e97ff7c0-2016-4de1-bea9-175c26381aa8)

  
  
## 🚀 How to Run and Examples

1. **Edit the `Makefile`** to update the path to IBM CPLEX and the path to the source files.  
2. **Open a terminal** in the folder where the `Makefile` is located.  
3. **Run the following commands:**  

```bash
chmod u+x makefile
./make
./FJS -i $instance -o $output -t $time-limit -a $learning-rate -s $seed
```
after seed parameter the other parameters are method-dependent.

The following methods do not need more parameters
  - Exact Models
	- Choose "MODELOMILP" in the sixth parameter for MILP Model
	- Choose "MODELOCP" in the sixth parameter for CP Model
	- Example Usage
    
		```
		./FJS -i DAFJS01.txt -o CP-Model-outout.csv -t 60 -a -0.3 -s 91287 -m MODELOCP
  		./FJS -i DAFJS01.txt -o MILP-Model-outout.csv -t 60 -a -0.3 -s 91287 -m MODELOMILP
		```
  
  - Tayebi Genetic Algorithm
  	- Choose "Tayebi" for the sixth parameter
   	- Example Usage
      
		```
		./FJS -i DAFJS01.txt -o Tayebi-outout.csv -t 60 -a -0.3 -s 91287 -m Tayebi
		```
    
  - Constructive heuristic
  	- Its value is given in the nineth parameters. Values can be ECT, SPT or Best for the best between ECT and SPT.
   	- Example Usage
      
    		```
      		./FJS -i DAFJS01.txt -o ECT-outout.csv -t 60 -a -0.3 -s 91287 -ls None -lse None -he ECT -mh None
      		```
      
  - Local Search
  	- Its value is given in the seventh (Local Search Neighborhood \in {Full, Reduced, CriticalReduced}) and eigthth (Local Search Strategy \in {Best Improvement, First Improvement}) parameters.
	- Example Usage
    
    		```
      		./FJS -i DAFJS01.txt -o LocalSearchReducedNeiBestImprovECT-outout.csv -t 60 -a -0.3 -s 91287 -ls Reduced -lse Best -he ECT -mh None
      		```
    
  - Metaheuristics:
   	For the metaheuristics we have the common parameters 10, 11 and 12, which stands for a tolerance, criticalOperations, and max number of iterations, respectively.
    	If tolerance is 0.1, the local search stops when it could improve at least 10% of the initial solution, if tolerance is 0.0, the local search will run until the end.
    	CriticalOperations is a bolean parameter, if it is set to 1, only critical operations will be considered in the local search (cropped neighborhood), else, all the operations will be considered in the local search (reduced neightbourhood).
    
	- ILS
    		It has parameters ellmin and ellmax
   
    		```
    		./FJS -i DAFJS01.txt -o SA-outout.csv -t 60 -a -0.3 -s 91287 -ls Reduced -lse Best -he Best -mh ILS -tol 0 -c 0 -itmax 1 -ellmin 2 -ellmax 4 &
    		```
    
	- Tabu Search 
		- tsize = size of Tabu List

  		```
    		./FJS -i DAFJS01.txt -o TS-outout.csv -t 60 -a -0.3 -s 91287 -ls Reduced -lse Best -he Best -mh TS -tol 0 -c 0 -itmax 0 -tsize 9 
    		```
    
	- Simulated Annealing
  		- It considers the parameters pertMin, pertMax, T0m, T0p, Tf, deltaMin, deltaMax.

  		- Example Usage:

      		```
    		./FJS -i DAFJS01.txt -o SA-outout.csv -t 60 -a -0.3 -s 91287 -ls Reduced -lse Best -he Best -mh SA -tol 0 -c 0 -itmax -1 --pertMin 3 --pertMax 3 --T0m 0.78 --T0p 0.79 --Tf 0.001 --deltaMin 0.82 --deltaMax 0.82
    		```
    
	- GRASP
  		- It has parameters alpha to control the randomness 
   
    		- Example Usage:
      
		```
    		./FJS -i DAFJS01.txt -o SA-outout.csv -t 60 -a -0.3 -s 91287 -ls Reduced -lse Best -he Best -mh GRASP -tol 0 -c 1 -itmax 1 -alphaGRASP 0.59 &
    		```
