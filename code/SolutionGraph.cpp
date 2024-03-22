#include "SolutionGraph.h"

#define NMAX 10000

bool VISITED[NMAX];
int PI[NMAX];
int indegree[NMAX];

SolutionGraph::SolutionGraph(Cenario& FJSFSLE)
{
	nNodes = FJSFSLE.nOperations + 2;
	nMachines = FJSFSLE.nMachines;
	AdjList = vector<list<SolutionNode>>(nNodes);

	for (int i = 0; i < nNodes; i++)
		indegree[i] = 0;

	for (auto& a : FJSFSLE.PrecedenceArcSet) {
		AdjList[a.source + 1].push_back(SolutionNode(a.tail + 1, -1, 0.0));
		indegree[a.tail + 1]++;
	}

	//dummy arcs
	for (int v = 1; v < nNodes - 1; v++) {
		if (AdjList[v].empty())
			AdjList[v].push_back(SolutionNode(nNodes - 1, -1, 0.0));

		if (indegree[v] == 0)
			AdjList[0].push_back(SolutionNode(v, -1, 0.0));
	}

	makespan = FJSFSLE.BIG_M;
	positions = vector< int>(nNodes, 0);
	criticalPath = vector< int>(nNodes, 0);
	R_in = vector<bool>(nNodes, 0);
	R_out = vector<bool>(nNodes, 0);
	uu = vector< int>(nNodes, 0);
	ff = vector< int>(nNodes, 0);
	st = vector<int>(nNodes, 0);
	weight = vector<int>(nNodes, 0);
	Machines = vector<vector< int>>(nMachines);
	isACriticalOp = vector<bool>(nNodes, false);
	isACriticalMachine = vector<bool>(FJSFSLE.nMachines, false);
	lastCriticalPosition = vector< int>(FJSFSLE.nMachines, 0);
	timeToGetSolution = 0.0;
	LocalSearchIterations = 0;
	nNeighbours = 0;
}

SolutionGraph::SolutionGraph()
{
	nNodes = 0;
	nMachines = 0;
	v = 0;
	AdjList = vector<list<SolutionNode>>();
	makespan = 999999;
	positions = vector< int>();
	criticalPath = vector< int>();
	R_in = vector<bool>();
	R_out = vector<bool>();
	uu = vector< int>();
	ff = vector< int>();
	st = vector<int>();
	weight = vector<int>();
	Machines = vector<vector< int>>();
	isACriticalOp = vector<bool>();
	isACriticalMachine = vector<bool>();
	lastCriticalPosition = vector< int>();
	timeToGetSolution = 0.0;
	LocalSearchIterations = 0;
	nNeighbours = 0;
}

void SolutionGraph::TopologicalSort() {
	for (int i = 0; i < nNodes; i++)
		indegree[i] = 0;
	for (int i = 0; i < this->nNodes; i++)
		for (auto& j : this->AdjList[i])
			indegree[j.id]++;

	queue<int> Queue;
	Queue.push(0);
	int iterator = 0;

	while (!Queue.empty()) {
		int i = Queue.front();
		Queue.pop();
		this->uu[iterator] = i; iterator++;

		for (auto& j : this->AdjList[i])
			if (--indegree[j.id] == 0)
				Queue.push(j.id);
	}
}

void SolutionGraph::TopologicalSortRecursive() {
	for (int i = 0; i < nNodes; i++)
		VISITED[i] = 0;
	int iterador = this->nNodes - 1;

	TopologicalSortUtil(0, this->uu, iterador);
}

void SolutionGraph::TopologicalSortUtil(const  int& i, vector< int>& ulist, int& iterador) {
	VISITED[i] = true;
	for (auto& j : AdjList[i]) {
		if (!VISITED[j.id]) {
			TopologicalSortUtil(j.id, ulist, iterador);
		}
	}
	ulist[iterador] = i;
	iterador--;
}

void SolutionGraph::TopologicalSortPlus(const  int& i, vector< int>& ulist, int& iterador)
{
	VISITED[i] = true;
	for (auto& j : AdjList[i]) {
		if (!VISITED[j.id])
			TopologicalSortPlus(j.id, this->uu, iterador);

		if (!R_in[i] && R_in[j.id])
			R_in[i] = true;
	}
	ulist[iterador] = i;
	iterador--;
}

void SolutionGraph::CriticalPath(Cenario& FJSSFLE, bool precisaReordenar) {
	FJSSFLE.nAvaliacoesMakespan++;
	fill(st.begin(), st.end(), -FJSSFLE.BIG_M);

	fill(isACriticalOp.begin(), isACriticalOp.end(), false);
	fill(isACriticalMachine.begin(), isACriticalMachine.end(), false);

	st[0] = 0;

	if (precisaReordenar)
		this->TopologicalSortRecursive();

	for (int ell = 0; ell < this->nNodes; ell++)
	{
		int i = this->uu[ell];

		for (auto& j : this->AdjList[i])
		{
			if (st[j.id] < st[i] + this->weight[i])
			{
				st[j.id] = st[i] + this->weight[i];
				PI[j.id] = i;
			}
		}
	}

	this->makespan = st[this->nNodes - 1];

	
	int i = 0;  
	if (this->nNodes > 0)
		i = PI[((int)this->nNodes - 1)];
	else {
		cerr << "ERRO" << endl;
		exit(0);
	}

	this->criticalPath.clear();

	for (int k = 0; k < FJSSFLE.nMachines; k++)
	{
		this->lastCriticalPosition[k] = 0;
	}

	int counterAntiLooping = 0;
	do
	{
		if ((int)this->ff[i] != -1) {
			if (this->lastCriticalPosition[this->ff[i]] == 0)
			{
				this->lastCriticalPosition[this->ff[i]] = this->positions[i];
			}
		}
		this->criticalPath.push_back(i);

		this->isACriticalOp[i] = true;

		if (this->weight[i] > 0) {
			this->isACriticalMachine[this->ff[i]] = true;
		}

		i = PI[i];

		if (++counterAntiLooping > nNodes)
		{
			cerr << "Infeasible solution - loop" << endl;
			exit(-1);
		}
	} while (i != 0);

	this->criticalPath = vector< int>(criticalPath.rbegin(), criticalPath.rend());
}

bool SolutionGraph::Reaches(const  int& i, const  int& j)
{
	for (int i = 0; i < nNodes; i++)
		VISITED[i] = 0;

	return ReachesUtil(i, j);
}

//se xx é alcançado por u
bool SolutionGraph::ReachesUtil(const  int& i, const  int& j)
{
	VISITED[i] = true;

	for (auto& v : AdjList[i])
	{
		if (v.id == j)
		{
			VISITED[j] = true;
			return true;
		}
		if (!VISITED[v.id])
		{
			ReachesUtil(v.id, j);
		}
	}
	return VISITED[j];
}

bool SolutionGraph::Feasible(Cenario& FJSSFLE) {
	for (int k = 0; k < FJSSFLE.nMachines; k++)
	{
		if ((int)this->Machines[k].size() == 0) continue;
		for (int i = 0; i < (int)this->Machines[k].size() - 1; i++)
		{
			int actualI = this->Machines[k][i];
			int actualInext = this->Machines[k][i + 1];
			if (this->st[actualI] + this->weight[actualI] > this->st[actualInext]) {
				return false;
			}
		}
	}

	for (int k = 0; k < nMachines; k++) {
		if ((int)Machines[k].size() > 1) {
			if ((int)this->Machines[k].size() == 0) continue;
			for (int i = 0; i < (int)Machines[k].size() - 1; i++) {
				if (AdjList[Machines[k][i]].back().id != Machines[k][i + 1]) {
					//cout << "\n>>>>>> Infeasible" << endl;
					return false;
				}
			}
		}
	}

	for (int i = 0; i < nNodes; i++) {
		int cont = 0;
		for (auto& adj : AdjList[i]) {
			if (adj.label >= 0) {
				cont++;
				if (cont != 1) {
					// cout << "\n>>>>>> Infeasible" << endl;
					return false;
				}
			}
		}
	}

	for (int i = 1; i < nNodes - 1; i++) {
		if (find(FJSSFLE.F[i - 1].begin(), FJSSFLE.F[i - 1].end(), ff[i]) == FJSSFLE.F[i - 1].end()) {
			//cout << "\n>>>>>> Infeasible" << endl;
			return false;
		}
	}

	//cout << "\n>>>>>> Solution graph is feasible" << endl;
	return true;
}

void SolutionGraph::DFSPlus(int i) {
	R_out[i] = true;
	for (auto& j : AdjList[i])
		if (R_out[j.id] == false)
			DFSPlus(j.id);
}

void SolutionGraph::RemoveOperation(Cenario& FJSSFLE, int v)
{
	this->v = v;
	weight[v] = 0;

	//remover arco de saída
	//remover o último arco pq será o único de máquina
	int suc = -1;
	if (this->AdjList[v].back().label != -1)
	{
		suc = this->AdjList[v].back().id;
		this->AdjList[v].pop_back();
	}

	//remover arco de entrada se houver
	if (positions[v] > 1) {
		int prec = this->Machines[this->ff[v]][this->positions[v] - 2];
		if (this->AdjList[prec].back().label != -1)
			AdjList[prec].pop_back();
		//religar operações na máquina com ausencia de v
		if (suc != -1)
			AdjList[prec].push_back(SolutionNode(suc, this->ff[v], this->weight[suc]));
	}

	/*
	* Atualizar posições, tempos de processamento e conjunto de operações na máquina
	*/
	int v_position = positions[v];
	positions[v] = 0;
	for (int i = v_position; i < (int)Machines[ff[v]].size(); i++)
	{
		int operation = Machines[ff[v]][i];
		weight[operation] = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[operation - 1][ff[v]], i, FJSSFLE.learningRate);
		positions[operation]--;
		Machines[ff[v]][i - 1] = Machines[ff[v]][i];
	}
	Machines[ff[v]].pop_back();
	ff[v] = -1;
	for (int i = 0; i < nNodes; i++)
		VISITED[i] = 0;

	fill(this->R_in.begin(), this->R_in.end(), false);
	R_in[this->v] = true;

	int iterador = this->nNodes - 1;

	TopologicalSortPlus(0, this->uu, iterador);

	fill(this->R_out.begin(), this->R_out.end(), false);

	DFSPlus(v);

	CriticalPath(FJSSFLE, false);
}

void SolutionGraph::InsertOperation(Cenario& FJSSFLE, int v, int gamma, int kappa) {
	ff[v] = kappa;
	positions[v] = gamma;
	Machines[kappa].push_back(v);

	for (int i = (int)Machines[kappa].size(); i-- > (gamma - 1); )
	{
		if (i == gamma - 1)
			Machines[kappa][i] = v;
		else {
			Machines[kappa][i] = Machines[kappa][i - 1];
			positions[Machines[kappa][i]] = i + 1;
		}

		int operation = Machines[kappa][i];
		weight[operation] = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[operation - 1][kappa], positions[operation], FJSSFLE.learningRate);
	}
	if (positions[v] > 1)
	{
		if (AdjList[Machines[kappa][positions[v] - 2]].back().label != -1)
			AdjList[Machines[kappa][positions[v] - 2]].pop_back();

		AdjList[Machines[kappa][positions[v] - 2]].push_back(SolutionNode(v, kappa, weight[v]));
	}
	if (positions[v] < (int)Machines[kappa].size()) {
		AdjList[v].push_back(SolutionNode(Machines[kappa][positions[v]], kappa, weight[Machines[kappa][positions[v]]]));
	}

	for (int i = 0; i < nNodes; i++)
		VISITED[i] = 0;
	fill(this->R_in.begin(), this->R_in.end(), false);
	R_in[this->v] = true;

	int iterador = this->nNodes - 1;
	TopologicalSortPlus(0, this->uu, iterador);

	CriticalPath(FJSSFLE, false);
}

void SolutionGraph::Perturbate(Cenario& FJSSFLE)
{
	v = rand() % (nNodes - 2) + 1;

	RemoveOperation(FJSSFLE, v);

	int randomPosition = rand() % ((int)FJSSFLE.F[v - 1].size());
	int kappa = FJSSFLE.F[v - 1][randomPosition];

	//generate random feasible positions
	int
		gamma_lb = 1,
		gamma_ub = (int)Machines[kappa].size() + 1;

	vector< int> feasiblePositions;

	for (int ell = 0; ell < (int)Machines[kappa].size(); ell++)
	{
		int operation = Machines[kappa][ell];
		if (R_in[operation])
			gamma_lb = ell + 2;
		if (R_out[operation] && gamma_ub == (int)Machines[kappa].size() + 1)
			gamma_ub = ell + 1;
	}

	int gamma = (rand() % (gamma_ub - gamma_lb + 1)) + gamma_lb;
	InsertOperation(FJSSFLE, v, gamma, kappa);
}

void SolutionGraph::ReducedPerturbate(Cenario& FJSSFLE, bool onlycriticaloperations)
{
	if (onlycriticaloperations)
	{
		v = rand() % ((int)this->criticalPath.size()) + 1;
		v = this->criticalPath[v];
	}
	else
	{
		v = rand() % (nNodes - 2) + 1;
	}

	int MakespanG = this->makespan;
	int k0 = this->ff[v];
	int r0 = this->positions[v];

	RemoveOperation(FJSSFLE, v);

	int randomPosition = rand() % ((int)FJSSFLE.F[v - 1].size());
	int kappa = FJSSFLE.F[v - 1][randomPosition];

	//generate random feasible positions
	int
		gamma_lb = 1,
		gamma_ub = (int)Machines[kappa].size() + 1;

	for (int ell = 0; ell < (int)this->Machines[kappa].size(); ell++)
	{
		int operation = this->Machines[kappa][ell];
		if (this->R_in[operation])
			gamma_lb = this->positions[operation] + 1;
		if (this->R_out[operation] && gamma_ub == (int)this->Machines[kappa].size() + 1)
			gamma_ub = this->positions[operation];
	}

	if (this->makespan >= MakespanG)
		gamma_ub = min(gamma_ub, this->lastCriticalPosition[kappa]);

	if (gamma_ub < gamma_lb) {
		InsertOperation(FJSSFLE, v, r0, k0);
		return;
	}

	int gamma = (rand() % (gamma_ub - gamma_lb + 1)) + gamma_lb;
	InsertOperation(FJSSFLE, v, gamma, kappa);
}

void SolutionGraph::printGanttChart(string instanceName, string learningRate, string method, vector<int>& jobs) {
	string
		saidaGantt = instanceName,
		saidaSolucao = instanceName;
	saidaGantt += ".solu";
	saidaSolucao += ".out";
	string figureName = "output\\" + instanceName;
	size_t pos = figureName.find(".txt");
	if (pos != string::npos)
		figureName.erase(pos, 4);
	pos = figureName.find("instances\\");
	if (pos != string::npos)
		figureName.erase(pos, 10);
	figureName += "_" + learningRate + "_" + method;

	ofstream outfile(saidaSolucao);

	outfile << " makespan (z) found: " << makespan;
	outfile << "===\n";

	for (int i = 0; i < nNodes - 2; i++) {
		outfile << "s_{" << i + 1 << "} = " << st[i + 1] << "\n";
	}

	for (int i = 0; i < nNodes - 2; i++) {
		outfile << "p_{" << i + 1 << "} = " << weight[i + 1] << "\n";
	}

	for (int i = 0; i < nNodes - 2; i++) {
		outfile << "ff_{" << i + 1 << "} = " << ff[i + 1] << "\n";
	}

	outfile.close();

	vector<double> timeAccumulatedInMachines(nMachines, 0.0);

	outfile.open(saidaGantt);
	outfile << "1,0,90,Type 0,-1" << endl;
	for (int i = 1; i <= nNodes - 2; i++) {
		outfile << ff[i] << "," << st[i] / 100 << "," << (st[i] + weight[i]) / 100 << ",Type " << jobs[i - 1] << "," << i << endl;
		//outfile << ff[i] << "," << st[i] << "," << (st[i] + weight[i]) << ",Type " << jobs[i - 1] << "," << i << endl;
	}
	outfile.close();

	string
		one = "python -c \"from imprimir_gantt import imprimir_gantt; imprimir_gantt('";

	one += saidaGantt;
	one += "','";
	one += figureName;
	one += ".png')\"";

	//int result = system("python -c 'print('oi')'");
	int result = system(one.c_str());
	cout << one << " " << result << endl;
}