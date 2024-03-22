#include "Cenario.h"

int Cenario::learningFunction(int p, int pos, double alpha)
{
	return static_cast<int>(round((double)p * pow(pos, alpha)));
}

Cenario::Cenario(const char* inputName, const char* outpuName, double learningRate, const double timeLimit, int SEED) {
	this->seed = SEED;
	srand(seed);

	nAvaliacoesMakespan = 0;
	this->inputName = inputName;
	this->outputName = outpuName;
	this->learningRate = learningRate;
	this->timeLimit = timeLimit;

	ifstream file(inputName);
	if (file.fail())
	{
		cerr << "Error: input failed!" << endl;
		exit(-1);
	}
	bool areClassicInstances = true;
	if (areClassicInstances) {
		nOperations = 0;
		nArcs = 0;
		file >> nJobs >> nMachines;
		Jobs = vector<vector< int>>(nJobs);
		ProcessingTimes = vector<vector<int>>(1000);
		OperationSet = vector<Operation>(1000, Operation());
		F = vector<vector< int>>(1000);
		I = vector<vector< int>>(nMachines);
		AdjBack = vector<list<int>>(1000);
		AdjFront = vector<list<int>>(1000);

		for (int i = 0; i < 1000; i++)
			ProcessingTimes[i] = vector<int>(nMachines, 0);

		for (int j = 0; j < nJobs; j++)
		{
			int nOpPerJob = 0;
			file >> nOpPerJob;
			Jobs[j] = vector< int>(nOpPerJob);
			for (int i = 0; i < nOpPerJob; i++) {
				Jobs[j][i] = nOperations++;
				int operation = Jobs[j][i];
				int nMachinesAux = 0;
				file >> nMachinesAux;
				OperationSet[operation].idJob = j;
				OperationSet[operation].idOperation = operation;
				OperationSet[operation].nMachines = nMachinesAux;
				OperationSet[operation].machineSet = vector< int>(nMachinesAux);
				OperationSet[operation].processingTimes = vector<int>(nMachinesAux);

				for (int k = 0; k < nMachinesAux; k++)
				{
					int machineAux, processingTimeAux;
					file >> machineAux >> processingTimeAux;
					OperationSet[operation].machineSet[k] = machineAux;
					OperationSet[operation].processingTimes[k] = 100 * processingTimeAux;
					ProcessingTimes[operation][machineAux] = 100 * processingTimeAux;
					I[machineAux].push_back(operation);
					F[operation].push_back(machineAux);
				}
				if (i > 0) {
					PrecedenceArcSet.push_back(Arc(Jobs[j][i - 1], Jobs[j][i], 0));
					nArcs++;
					AdjBack[Jobs[j][i]].push_back(Jobs[j][i - 1]);
					AdjFront[Jobs[j][i - 1]].push_back(Jobs[j][i]);
				}
				std::sort(OperationSet[operation].machineSet.begin(), OperationSet[operation].machineSet.end());
				std::sort(F[operation].begin(), F[operation].end());

				//ordenar as m�quinas pelo indice
				vector< int> v = OperationSet[operation].machineSet;
				vector<size_t> idx(nMachinesAux);
				iota(idx.begin(), idx.end(), 0);
				stable_sort(idx.begin(), idx.end(),
					[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
				v = OperationSet[operation].machineSet;
				vector<int> vv = OperationSet[operation].processingTimes;

				for (int k = 0; k < nMachinesAux; k++)
				{
					OperationSet[operation].machineSet[k] = v[idx[k]];
					OperationSet[operation].processingTimes[k] = ProcessingTimes[operation][OperationSet[operation].machineSet[k]];
				}
			}
		}

		OperationSet.resize(nOperations);
		ProcessingTimes.resize(nOperations);
		F.resize(nOperations);
		AdjBack.resize(nOperations);
		AdjFront.resize(nOperations);
	}
	else {
		//se s�o inst�ncias do ernesto descomenta essa parte
		int aux1, aux2;
		file >> aux1 >> aux2;
		//---------------------------------------------------

		file >> nOperations >> nArcs >> nMachines;

		PrecedenceArcSet = vector<Arc>(nArcs);
		OperationSet = vector<Operation>(nOperations);
		I = vector<vector< int>>(nMachines);
		F = vector<vector< int>>(nOperations);
		ProcessingTimes = vector<vector<int>>(nOperations);

		for (int i = 0; i < nOperations; i++)
			ProcessingTimes[i] = vector<int>(nMachines);

		AdjBack = vector<list<int>>(nOperations);
		AdjFront = vector<list<int>>(nOperations);

		//ler arcos de precedencia
		for (int i = 0; i < nArcs; i++) {
			int source, tail;
			file >> source >> tail;
			PrecedenceArcSet[i] = Arc(source, tail, 0);

			AdjBack[tail].push_back(source);
			AdjFront[source].push_back(tail);
		}

		for (int i = 0; i < nOperations; i++) {
			int nMachinesAux;
			file >> nMachinesAux;
			OperationSet[i].nMachines = nMachinesAux;
			OperationSet[i].machineSet = vector< int>(nMachinesAux);
			OperationSet[i].processingTimes = vector<int>(nMachinesAux);

			for (int k = 0; k < nMachinesAux; k++) {
				int machineAux, processingTimeAux;
				file >> machineAux >> processingTimeAux;
				OperationSet[i].machineSet[k] = machineAux;
				OperationSet[i].processingTimes[k] = 100 * processingTimeAux;
				ProcessingTimes[i][machineAux] = 100 * processingTimeAux;
				I[machineAux].push_back(i);
				F[i].push_back(machineAux);
			}
			std::sort(OperationSet[i].machineSet.begin(), OperationSet[i].machineSet.end());
			std::sort(F[i].begin(), F[i].end());

			//ordenar as m�quinas pelo indice
			vector< int> v = OperationSet[i].machineSet;
			vector<size_t> idx(nMachinesAux);
			iota(idx.begin(), idx.end(), 0);
			stable_sort(idx.begin(), idx.end(),
				[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
			v = OperationSet[i].machineSet;
			vector<int> vv = OperationSet[i].processingTimes;

			for (int k = 0; k < nMachinesAux; k++)
			{
				OperationSet[i].machineSet[k] = v[idx[k]];
				OperationSet[i].processingTimes[k] = ProcessingTimes[i][OperationSet[i].machineSet[k]];
			}
		}
		vector<bool> isVisited = vector<bool>(nOperations, false);
		vector<list<int>> AuxAdjFront = AdjFront;
		for (int i = 0; i < nOperations; i++)
		{
			for (auto& j : AdjFront[i]) {
				AuxAdjFront[j].push_back(i);
			}
		}

		int jobCount = 0;
		for (int i = 0; i < nOperations; i++) {
			if (!isVisited[i]) {
				isVisited[i] = true;
				OperationSet[i].idJob = jobCount;
				jobCount++;
				CountJobs(isVisited, AuxAdjFront, i, OperationSet, jobCount);
			}
		}
		nJobs = jobCount;
		Jobs = vector<vector< int>>(nJobs);
		for (int i = 0; i < nOperations; i++) {
			OperationSet[i].idJob--;
			Jobs[OperationSet[i].idJob].push_back(i);
		}
	}
	file.close();

	// compute value of infinite
	BIG_M = 1;
	for (int v = 0; v < nOperations; ++v) {
		int largest = 0;
		for (int k = 0; k < OperationSet[v].nMachines; ++k)
			if (largest < OperationSet[v].processingTimes[k])
				largest = OperationSet[v].processingTimes[k];
		BIG_M += largest;
	}
	BIG_M = INT_MAX;

	int min = 10000, max = 0, sum = 0;
	for (auto& i : Jobs)
	{
		if ((int)i.size() < min)
			min = (int)i.size();
		if ((int)i.size() > max)
			max = (int)i.size();
	}
	vector<int> Machines(nMachines, 0);
	for (auto& i : this->OperationSet)
	{
		for (auto& j : i.machineSet)
		{
			Machines[j]++;
		}
	}
	for (auto& i : Machines)
		sum += i;

	double omega_1 = 0.0;

	for (auto& job : Jobs)
	{
		double  a_kappa = 0,
			a_kappa_min = (double)job.size() - 1,
			a_kappa_max = (double)job.size() * (job.size() - 1) / 2.0;
		for (auto& i : job) {
			queue<int> Queue;
			vector<bool> visited(nOperations, false);
			Queue.push(i);
			while (!Queue.empty()) {
				for (auto& adj : AdjFront[Queue.front()]) {
					Queue.push(adj);

					if (!visited[adj])
						a_kappa++;

					visited[adj] = true;
				}
				Queue.pop();
			}
		}
		omega_1 += 1 - (a_kappa - a_kappa_min) / (a_kappa_max - a_kappa_min);
	}
	omega_1 /= (double)Jobs.size();

	double omega_2 = 0.0;

	for (auto& machine : I)
	{
		omega_2 += machine.size();
	}
	omega_2 = (double)(omega_2 - nOperations) / (double)(nOperations * nMachines - nOperations);

	ofstream saida("caracteristicas.csv", fstream::app);

	saida << this->inputName << "," << this->nMachines << "," << this->nOperations << "," << this->nJobs << "," << this->nArcs << "," << sum << "," << omega_1 << "," << omega_2 << endl;

	saida.close();
}