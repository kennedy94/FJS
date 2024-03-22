#include "FJS.h"

int learningFunction(int p, int pos, double alpha)
{
	return (int)round((double)p * pow(pos, alpha));
}

//se xx é alcançado por u
//se xx é alcançado por u
bool alcancado(vector<bool>& visited, vector<vector<int>>& SJ, int u, int xx) {
	visited[u] = true;
	for (auto& v : SJ[u]) {
		if (v == xx) {
			visited[xx] = true;
			return true;
		}

		if (!visited[v]) {
			alcancado(visited, SJ, v, xx);
		}
	}
	return visited[xx];
}

bool alcanca(vector<bool>& visited, vector<vector<int>>& SJ, int u, int xx) {
	for (int i = 0; i < (int)visited.size(); i++) {
		visited[i] = false;
	}
	return alcancado(visited, SJ, u, xx);
}

void DFS_count_job(vector<bool>& visited, vector<list<int>> AdjList, int vertice, vector<operation>& op, int job) {
	op[vertice].job = job;
	visited[vertice] = true;
	for (auto& adj : AdjList[vertice])
		if (!visited[adj])
			DFS_count_job(visited, AdjList, adj, op, job);
}

FJS::FJS(const char* instance, const char* saida, double BUDGET, double alpha) {
	this->instanceName = instance;
	this->outputName = saida;
	this->alpha = alpha;
	this->TimeLimit = BUDGET;
	this->SEED = (int)time(NULL);
	srand(this->SEED);
	ifstream FILE(instance);

	//se são instâncias do ernesto descomenta essa parte
	//int inu1, inu2;
	//FILE >> inu1 >> inu2;
	//---------------------------------------------------
	bool areClassicInstances = true;
	if (areClassicInstances) {
		nOp = 0;
		nArcs = 0;
		FILE >> nJobs >> nMac;
		Jobs = vector<vector< int>>(nJobs);
		processingTimesMatrix = vector<vector<int>>(1000);
		operationsVector = vector<operation>(1000);
		//F = vector<vector< int>>(1000);
		I = vector<vector< int>>(nMac);
		AdjBack = vector<vector<int>>(1000);
		AdjFront = vector<vector<int>>(1000);

		for (int i = 0; i < 1000; i++)
			processingTimesMatrix[i] = vector<int>(nMac, -1);

		for (int j = 0; j < nJobs; j++)
		{
			int nOpPerJob = 0;
			FILE >> nOpPerJob;
			Jobs[j] = vector< int>(nOpPerJob);
			for (int i = 0; i < nOpPerJob; i++) {
				Jobs[j][i] = nOp++;
				int operation = Jobs[j][i];
				int nMachinesAux = 0;
				FILE >> nMachinesAux;
				operationsVector[operation].job = j;
				operationsVector[operation].Fsize = nMachinesAux;
				operationsVector[operation].F = vector< int>(nMachinesAux);
				operationsVector[operation].P = vector<int>(nMachinesAux);

				for (int k = 0; k < nMachinesAux; k++)
				{
					int machineAux, processingTimeAux;
					FILE >> machineAux >> processingTimeAux;
					operationsVector[operation].F[k] = machineAux;
					operationsVector[operation].P[k] = 100 * processingTimeAux;
					processingTimesMatrix[operation][machineAux] = 100 * processingTimeAux;
					I[machineAux].push_back(operation);
					//F[operation].push_back(machineAux);
				}
				if (i > 0) {
					precedenceArcsVector.push_back(arc(Jobs[j][i - 1], Jobs[j][i]));
					nArcs++;
					AdjBack[Jobs[j][i]].push_back(Jobs[j][i - 1]);
					AdjFront[Jobs[j][i - 1]].push_back(Jobs[j][i]);
				}
				std::sort(operationsVector[operation].F.begin(), operationsVector[operation].F.end());
				//std::sort(F[operation].begin(), F[operation].end());

				//ordenar as m�quinas pelo indice
				vector< int> v = operationsVector[operation].F;
				vector<size_t> idx(nMachinesAux);
				iota(idx.begin(), idx.end(), 0);
				stable_sort(idx.begin(), idx.end(),
					[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
				v = operationsVector[operation].F;
				vector<int> vv = operationsVector[operation].P;

				for (int k = 0; k < nMachinesAux; k++)
				{
					operationsVector[operation].F[k] = v[idx[k]];
					operationsVector[operation].P[k] = processingTimesMatrix[operation][operationsVector[operation].F[k]];
				}
			}
		}

		operationsVector.resize(nOp);
		processingTimesMatrix.resize(nOp);
		//F.resize(nOperations);
		AdjBack.resize(nOp);
		AdjFront.resize(nOp);
	}
	else{
		FILE >> nOp >> nArcs >> nMac;

		precedenceArcsVector = vector<arc>(nArcs);
		operationsVector = vector<operation>(nOp);
		I = vector<vector<int>>(nMac);
		processingTimesMatrix = vector<vector<int>>(nOp);
		for (int i = 0; i < nOp; i++) {
			processingTimesMatrix[i] = vector<int>(nMac, -1);
		}

		AdjBack = vector<vector<int>>(nOp);
		AdjFront = vector<vector<int>>(nOp);

		//ler arcos de precedencia
		for (int i = 0; i < nArcs; i++) {
			FILE >> precedenceArcsVector[i].i >> precedenceArcsVector[i].j;
			AdjBack[precedenceArcsVector[i].j].push_back(precedenceArcsVector[i].i);
			AdjFront[precedenceArcsVector[i].i].push_back(precedenceArcsVector[i].j);
		}

		for (int i = 0; i < nOp; i++) {
			FILE >> operationsVector[i].Fsize;

			operationsVector[i].F = vector<int>(operationsVector[i].Fsize);
			operationsVector[i].P = vector<int>(operationsVector[i].Fsize);

			for (int k = 0; k < operationsVector[i].Fsize; k++) {
				FILE >> operationsVector[i].F[k] >> operationsVector[i].P[k];
				operationsVector[i].P[k] *= 100;
				processingTimesMatrix[i][operationsVector[i].F[k]] = operationsVector[i].P[k];

				I[operationsVector[i].F[k]].push_back(i);
			}

			vector<int> v = operationsVector[i].F;

			vector<size_t> idx(operationsVector[i].Fsize);
			iota(idx.begin(), idx.end(), 0);
			stable_sort(idx.begin(), idx.end(),
				[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
			v = operationsVector[i].F;
			for (int k = 0; k < operationsVector[i].Fsize; k++)
				operationsVector[i].F[k] = v[idx[k]];

			vector<int> vv = operationsVector[i].P;
			for (int k = 0; k < operationsVector[i].Fsize; k++)
				operationsVector[i].P[k] = vv[idx[k]];
		}
	}
	FILE.close();

	vector<bool> isVisited(nOp, false);

	relacao_precedencia = vector<vector<int>>(nOp);
	for (int i = 0; i < nOp; i++) {
		relacao_precedencia[i] = vector<int>(nOp);
		for (int j = 0; j < nOp; j++) {
			if (i == j) {
				relacao_precedencia[i][j] = 0;
			}
			else {
				if (alcanca(isVisited, AdjFront, i, j)) {
					relacao_precedencia[i][j] = -1;
				}
				else {
					if (alcanca(isVisited, AdjBack, i, j)) {
						relacao_precedencia[i][j] = 1;
					}
					else {
						relacao_precedencia[i][j] = 0;
					}
				}
			}
		}
	}

	{
		vector<list<int>> AdjList(nOp);

		for (auto& i_arc : precedenceArcsVector) {
			AdjList[i_arc.i].push_back(i_arc.j);
			AdjList[i_arc.j].push_back(i_arc.i);
		}

		isVisited = vector<bool>(nOp, false);
		int jobCount = 0;
		for (int i = 0; i < nOp; i++) {
			if (!isVisited[i]) {
				isVisited[i] = true;
				operationsVector[i].job = jobCount;
				jobCount++;
				DFS_count_job(isVisited, AdjList, i, operationsVector, operationsVector[i].job);
			}
		}
		nJobs = jobCount;
		Jobs = vector<vector<int>>(nJobs);
		for (int i = 0; i < nOp; i++)
			Jobs[operationsVector[i].job].push_back(i);
	}
	// compute value of infinite
	infinite = 1;
	for (int v = 0; v < nOp; ++v) {
		int largest = 0;
		for (int k = 0; k < nMac; ++k)
			if (processingTimesMatrix[v][k] != -1.0 && largest < processingTimesMatrix[v][k])
				largest = processingTimesMatrix[v][k];
		infinite += largest;
	}

	for (int v = 0; v < nOp; ++v)
		for (int k = 0; k < nMac; ++k)
			if (processingTimesMatrix[v][k] == -1)
				processingTimesMatrix[v][k] = infinite;
}