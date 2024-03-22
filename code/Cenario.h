#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <fstream>
#include <numeric>
#include <algorithm>
#include "limits.h"
#include <math.h>
#include "GraphAlgorithms.h"
#include "Operation.h"
#include "Arc.h"

using namespace std;

class Cenario
{
public:
	int learningFunction(
		int p,   //normal processing time
		int pos,    //position in the machine
		double alpha //learning rate
	);

	int nOperations,
		nArcs,
		nMachines,
		nJobs;
	double timeLimit;
	int nAvaliacoesMakespan;

	vector<Operation> OperationSet;
	vector<Arc> PrecedenceArcSet;

	vector<vector< int>> I;	/*set of operations that
							can be produced by machines (nmachine x nop)
							*/
	vector<vector< int>> F;	/*set of machines that
							produce operations (nop x nmachine)
							*/
	vector<vector<int>> ProcessingTimes;	//processing times (nop x nmachine), dense matrix

	vector<list<int>> AdjBack, AdjFront; /*set of operations that directly precedes and follows,
								*respectively v
								*/

	vector<vector< int>> Jobs; /*
							   *Operações de cada job em ordem de indice
							   */

	string inputName, outputName;	//instance and output filenames

	double learningRate;	//time budget and learning rate
	int BIG_M;
	int seed;

	Cenario(const char* inputName, const char* outpuName, const double learningRate, const double timeLimit, int SEED);

	Cenario();
};
