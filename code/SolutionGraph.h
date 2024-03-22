#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <chrono>
#include <random>
#include <algorithm>
#include "SolutionNode.h"
#include "Cenario.h"

using namespace std;

class SolutionGraph
{
private:
	void TopologicalSort();
	void TopologicalSortRecursive();
	void TopologicalSortUtil(const  int& i, vector< int>& ulist, int& iterador);
	void TopologicalSortPlus(const  int& i, vector< int>& ulist, int& iterador);
public:

	int nNodes, nMachines, v, makespan;
	int LocalSearchIterations;
	int nNeighbours;
	double timeToGetSolution;

	vector<list<SolutionNode>> AdjList;

	vector<vector< int>> Machines;

	vector< int> positions,
		criticalPath,
		uu,
		ff,
		lastCriticalPosition;

	vector<bool>
		R_in,
		R_out;

	vector<int>
		st,
		weight;

	vector<bool>
		isACriticalOp,
		isACriticalMachine;

	SolutionGraph();

	SolutionGraph(Cenario& FJSFSLE);

	void CriticalPath(Cenario& FJSSFLE, bool precisaReordenar = true);
	bool Reaches(const  int& i, const  int& j);
	bool ReachesUtil(const  int& i, const  int& j);
	bool Feasible(Cenario& FJSSFLE);
	void DFSPlus(int i);
	void RemoveOperation(Cenario& FJSSFLE, int  v);
	void InsertOperation(Cenario& FJSSFLE, int v, int gamma, int kappa);

	void Perturbate(Cenario& FJSSFLE);

	void ReducedPerturbate(Cenario& FJSSFLE, bool onlycriticaloperations);

	SolutionGraph(const SolutionGraph& G)
	{
		nNodes = G.nNodes;
		nMachines = G.nMachines;
		AdjList = G.AdjList;
		v = G.v;
		makespan = G.makespan;
		positions = G.positions;
		criticalPath = G.criticalPath;
		R_in = G.R_in;
		R_out = G.R_out;
		uu = G.uu;
		ff = G.ff;
		st = G.st;
		weight = G.weight;
		Machines = G.Machines;
		isACriticalOp = G.isACriticalOp;
		isACriticalMachine = G.isACriticalMachine;
		lastCriticalPosition = G.lastCriticalPosition;
		timeToGetSolution = G.timeToGetSolution;
		LocalSearchIterations = G.LocalSearchIterations;
		nNeighbours = G.nNeighbours;
	}

	SolutionGraph& operator=(const SolutionGraph& G) {
		nNodes = G.nNodes;
		nMachines = G.nMachines;
		AdjList = G.AdjList;
		v = G.v;
		makespan = G.makespan;
		positions = G.positions;
		criticalPath = G.criticalPath;
		R_in = G.R_in;
		R_out = G.R_out;
		uu = G.uu;
		ff = G.ff;
		st = G.st;
		weight = G.weight;
		Machines = G.Machines;
		isACriticalOp = G.isACriticalOp;
		isACriticalMachine = G.isACriticalMachine;
		lastCriticalPosition = G.lastCriticalPosition;
		timeToGetSolution = G.timeToGetSolution;
		LocalSearchIterations = G.LocalSearchIterations;
		nNeighbours = G.nNeighbours;
		return *this;
	}

	void printGanttChart(string instanceName, string learningRate, string method, vector<int>& jobs);
};
