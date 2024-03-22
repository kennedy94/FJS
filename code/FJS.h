#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <fstream>
#include <numeric>      // std::iota
#include "SolutionGraph.h"

using namespace std;

int learningFunction(int p, int pos, double alpha);

struct operation {
	int Fsize = 0;
	int job = 0;
	vector<int> F;
	vector<int> P;
};

struct arc {
	int i, j;

	arc() {
		i = -1;
		j = -1;
	}
	arc(int i, int j) {
		this->i = i;
		this->j = j;
	}
};

/*
Class that represents all data from a given instance
*/
class FJS
{
public:
	//it receives file name for instance and output
	//time budget
	//learning rate
	FJS(const char* instance, const char* saida, double BUDGET, double alpha);

	int nOp,	//#operations
		nArcs,	//#precedence arcs
		nMac,	//#machines
		nJobs;	//number of jobs
	vector<vector<int>> I;	/*set of operations that
							can be produced by machines (nmachine x nop)
							*/
	vector<operation> operationsVector;	//vector of operations

	vector<arc> precedenceArcsVector;		//precedence arcs

	vector<vector<int>> processingTimesMatrix;	//processing times (nop x nmachine), dense matrix

	vector<vector<int>> AdjBack, AdjFront; /*set of operations that directly precedes and follows,
								*respectively v
								*/
	vector<vector<int>>
		relacao_precedencia; /*-1 se i vem antes de j no mesmo job,
							 * 0 se i não tem relação de precedência com j
							 * 1 se i vem depois de j no mesmo job
							 */

	vector<vector<int>> Jobs; /*
							   *Operações de cada job em ordem de indice
							   */

	string instanceName, outputName;	//instance and output filenames
	int SEED;

	double TimeLimit, alpha;	//time budget and learning rate
	int infinite;
	int nAvaliacoesMakespan = 0;
	int LocalSearchID, LSStrategyID;
};
