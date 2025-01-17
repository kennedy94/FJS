#include "Heuristics.h"

#define NMAX 10000

bool isCriticalOpG[NMAX];
int GTopologicalSorting[NMAX];

struct Moviment
{
	int i, k, ell;

	bool operator==(Moviment x) {
		return this->i == x.i && this->k == x.k;// && this->ell == x.ell;
	}
	Moviment() {
		this->i = -1;
		this->k = -1;
		this->ell = -1;
	};

	Moviment(int i, int k, int ell) {
		this->i = i;
		this->k = k;
		this->ell = ell;
	}

	Moviment(int i, int k) {
		this->i = i;
		this->k = k;
		this->ell = -1;
	}
};

SolutionGraph SPT(Cenario& FJSSFLE) {
	SolutionGraph G = SolutionGraph(FJSSFLE);

	vector<int>
		r_op(G.nNodes, FJSSFLE.BIG_M),
		r_mac(FJSSFLE.nMachines, 0),
		c(G.nNodes, FJSSFLE.BIG_M);
	r_op[0] = 0;
	c[0] = 0;

	vector< int>
		g(FJSSFLE.nMachines, 1);

	vector<bool> isScheduled(G.nNodes, false);
	isScheduled[0] = true;

	set<int> notScheduled;
	for (int i = 1; i < G.nNodes - 1; i++)
		notScheduled.insert(i);

	while (!notScheduled.empty())
	{
		list<int> readyOperations;

		int r_min = FJSSFLE.BIG_M;
		for (auto& v : notScheduled)
		{
			bool isReady = true;
			int readyTime = 0;
			for (auto& i : FJSSFLE.AdjBack[(v - 1)]) {
				if (!isScheduled[(i + 1)]) {
					isReady = false;
					break;
				}
				if (c[(i + 1)] > readyTime)
					readyTime = c[(i + 1)];
			}
			if (isReady) {
				r_op[v] = readyTime;
				readyOperations.push_back(v);

				for (auto& k : FJSSFLE.F[(v - 1)])
				{
					if (r_min > max(r_op[v], r_mac[k]))
						r_min = max(r_op[v], r_mac[k]);
				}
			}
		}

		int v_hat = -1, k_hat = -1;
		int  minProcessingtime = FJSSFLE.BIG_M;
		for (auto& v : readyOperations)
		{
			for (auto& k : FJSSFLE.F[v - 1])
			{
				if (r_min == max(r_op[v], r_mac[k]))
				{
					int realProcessingTime = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[v - 1][k], g[k], FJSSFLE.learningRate);

					if (realProcessingTime < minProcessingtime)
					{
						v_hat = v;
						k_hat = k;
						minProcessingtime = realProcessingTime;
					}
				}
			}
		}
		G.weight[v_hat] = minProcessingtime;
		G.ff[v_hat] = k_hat;
		c[v_hat] = r_min + minProcessingtime;
		r_mac[k_hat] = c[v_hat];
		G.positions[v_hat] = g[k_hat];
		G.st[v_hat] = r_min;
		g[k_hat]++;

		if (!G.Machines[k_hat].empty())
		{
			G.AdjList[G.Machines[k_hat].back()].push_back(SolutionNode(v_hat, k_hat, minProcessingtime));
		}
		G.Machines[k_hat].push_back(v_hat);
		notScheduled.erase(v_hat);
		isScheduled[v_hat] = true;
		//cout << "(" << v_hat << "," << k_hat << ")" << endl;
	}
	G.CriticalPath(FJSSFLE);
	return G;
}

SolutionGraph ECT(Cenario& FJSSFLE) {
	SolutionGraph G = SolutionGraph(FJSSFLE);

	vector<int>
		r_op(G.nNodes, FJSSFLE.BIG_M),
		r_mac(FJSSFLE.nMachines, 0),
		c(G.nNodes, FJSSFLE.BIG_M);
	r_op[0] = 0;
	c[0] = 0;

	vector< int>
		g(FJSSFLE.nMachines, 1);

	vector<bool> isScheduled(G.nNodes, false);
	isScheduled[0] = true;

	set<int> notScheduled;
	for (int i = 1; i < G.nNodes - 1; i++)
		notScheduled.insert(i);

	while (!notScheduled.empty())
	{
		int v_hat = -1, k_hat = -1;
		int
			minCompletionTime = FJSSFLE.BIG_M,
			chosenProcessingTime = FJSSFLE.BIG_M,
			chosenStartingTime = FJSSFLE.BIG_M;

		for (auto& v : notScheduled)
		{
			bool isReady = true;
			int readyTime = 0;
			for (auto& i : FJSSFLE.AdjBack[v - 1]) {
				if (!isScheduled[i + 1]) {
					isReady = false;
					break;
				}
				if (c[i + 1] > readyTime)
					readyTime = c[i + 1];
			}

			if (!isReady) continue;

			r_op[v] = readyTime;

			for (auto& k : FJSSFLE.F[v - 1])
			{
				int realProcessingTime = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[v - 1][k], g[k], FJSSFLE.learningRate);
				int startingTime = max(r_op[v], r_mac[k]);
				if (startingTime + realProcessingTime < minCompletionTime)
				{
					v_hat = v;
					k_hat = k;
					minCompletionTime = startingTime + realProcessingTime;
					chosenProcessingTime = realProcessingTime;
					chosenStartingTime = startingTime;
				}
			}
		}

		G.weight[v_hat] = chosenProcessingTime;
		G.ff[v_hat] = k_hat;
		c[v_hat] = minCompletionTime;
		r_mac[k_hat] = c[v_hat];
		G.positions[v_hat] = g[k_hat];
		G.st[v_hat] = chosenStartingTime;
		g[k_hat]++;

		if (!G.Machines[k_hat].empty())
		{
			G.AdjList[G.Machines[k_hat].back()].push_back(SolutionNode(v_hat, k_hat, chosenProcessingTime));
		}
		G.Machines[k_hat].push_back(v_hat);
		notScheduled.erase(v_hat);
		isScheduled[v_hat] = true;
		//cout << "(" << v_hat << "," << k_hat << ")" << endl;
	}
	G.CriticalPath(FJSSFLE);
	return G;
}

SolutionGraph LocalSearch(Cenario& FJSSFLE, SolutionGraph& G, string localSearchStrategy, double tolerancy, std::chrono::high_resolution_clock::time_point tic, bool onlyCriticalOperations) {
	SolutionGraph Gbest(G);
	double currentMakespan = Gbest.makespan;
#ifdef _DEBUG
	cout << "Local search:" << currentMakespan << endl;
#endif
	int iterations = 0,
		vizinhos = 0;
	do
	{
		iterations++;
		currentMakespan = Gbest.makespan;
		Gbest = ReducedLocalSearch(FJSSFLE, Gbest, localSearchStrategy, tic, onlyCriticalOperations);
		vizinhos += Gbest.nNeighbours;
#ifdef _DEBUG
		cout << "\t" << Gbest.makespan << endl;
#endif
	} while (-(Gbest.makespan - currentMakespan) / currentMakespan > tolerancy);
#ifdef _DEBUG
	cout << "fim" << endl;
#endif // _DEBUG
	Gbest.LocalSearchIterations = iterations;
	Gbest.nNeighbours = vizinhos;
	return Gbest;
}

SolutionGraph ReducedLocalSearch(Cenario& FJSSFLE, SolutionGraph& G, string localSearchStrategy, std::chrono::high_resolution_clock::time_point tic, bool onlyCriticalOperations)
{
	SolutionGraph Gbest = G;
	int nNodes = G.nNodes;
	for (int i = 0; i < nNodes; i++) {
		isCriticalOpG[i] = G.isACriticalOp[i];
		GTopologicalSorting[i] = G.uu[i];
	}
	int vizinhos = 0;
	int MakespanG = G.makespan;

	for (int i = 0; i < nNodes; i++)
	{
		int v = GTopologicalSorting[i];
		if (onlyCriticalOperations && !isCriticalOpG[v])
			continue;

		if (v == 0 || v == nNodes - 1)
			continue;

		int ellAntigo = G.positions[v],
			kAntigo = G.ff[v];
		G.RemoveOperation(FJSSFLE, v);

		for (auto& k : FJSSFLE.F[v - 1])
		{
			int
				gamma_lb = 1,
				gamma_ub = (int)G.Machines[k].size() + 1;

			for (int ell = 0; ell < (int)G.Machines[k].size(); ell++)
			{
				int operation = G.Machines[k][ell];
				if (G.R_in[operation])
					gamma_lb = G.positions[operation] + 1;
				if (G.R_out[operation] && gamma_ub == (int)G.Machines[k].size() + 1)
					gamma_ub = G.positions[operation];
			}

			if (G.makespan >= MakespanG)
				gamma_ub = min(gamma_ub, G.lastCriticalPosition[k]);

			for (int ell = gamma_lb; ell <= gamma_ub; ell++)
			{
				G.InsertOperation(FJSSFLE, v, ell, k);
				vizinhos++;
				if (Gbest.makespan > G.makespan) {
					Gbest = G;
					auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
					Gbest.timeToGetSolution = toc;
					if (localSearchStrategy == "First") {
						Gbest.nNeighbours = vizinhos;
						return Gbest;
					}
				}
				G.RemoveOperation(FJSSFLE, v);
			}
		}
		G.InsertOperation(FJSSFLE, v, ellAntigo, kAntigo);
	}
	Gbest.nNeighbours = vizinhos;
	return Gbest;
}

SolutionGraph ILS(Cenario& FJSSFLE, string localSearchStrategy, int pertMin, int pertMax, double tolerancy, bool tuning, bool onlyCriticalOperations) {
	auto tic = std::chrono::high_resolution_clock::now();

	SolutionGraph G, Gbest, GSPT, GECT, Gbest_;

	GSPT = SPT(FJSSFLE);
	GECT = ECT(FJSSFLE);

	G = (GSPT.makespan < GECT.makespan) ? GSPT : GECT;
	Gbest = G;

	string
		saidaDetalhada = "detailed" + FJSSFLE.outputName;

	size_t pos = saidaDetalhada.find(".txt");
	if (pos != string::npos)
		saidaDetalhada.erase(pos, 4);
	pos = saidaDetalhada.find("instances\\");
	if (pos != string::npos)
		saidaDetalhada.erase(pos, 10);

	int itSemMelhora = 0;

	{
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "ILS," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << 0 << "," << FJSSFLE.seed << endl;
		outfile.close();
	}

	int i = 1;
	double ultima_melhora = 0.0;

	while (true)
	{
		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
		if (toc >= FJSSFLE.timeLimit)
			break;
		if (tuning && (toc - ultima_melhora) > 5.0)
			break;

		G = LocalSearch(FJSSFLE, G, localSearchStrategy, tolerancy, tic, onlyCriticalOperations);

		if (G.makespan < Gbest.makespan) {
			Gbest = G;
			itSemMelhora = 0;
			auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
			Gbest.timeToGetSolution = toc;

			ofstream outfile(saidaDetalhada, fstream::app);
			outfile << "ILS," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << toc << "," << FJSSFLE.seed << endl;
			outfile.close();
			ultima_melhora = toc;
		}
		else {
			itSemMelhora++;
		}

		toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;

		int pertMaxCurrent = (rand() % (pertMax - pertMin + 1)) + pertMin;
		for (int j = 0; j < pertMaxCurrent; j++) {
			G.Perturbate(FJSSFLE);
		}
		i++;
	}

	{
		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "ILS," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << toc << "," << FJSSFLE.seed << endl;
		outfile.close();
	}

	return Gbest;
}

SolutionGraph SA(Cenario& FJSSFLE, string localSearch, string localSearchStrategy, int itMin, int itMax, double tempM, double tempP, double tempF, double deltaMin, double deltaMax, double tolerancy, string functionType, bool tuning)
{
	auto tic = std::chrono::high_resolution_clock::now();

	SolutionGraph G, Gbest, GSPT, GECT, GlocalSearch, Gline;

	GSPT = SPT(FJSSFLE);
	GECT = ECT(FJSSFLE);

	G = (GSPT.makespan < GECT.makespan) ? GSPT : GECT;
	Gbest = G;
	double initialMakespan = G.makespan;

	string
		saidaDetalhada = "detailed" + FJSSFLE.outputName;

	double temp = -tempM / log(tempP);
	double temp0 = temp;

	int itSemMelhora = 0;
	double ultima_melhora = 0.0;

	{
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "SA," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << 0 << "," << FJSSFLE.seed << endl;
		outfile.close();
	}

	while (true)
	{
		double tempFexp = tempF;
		double temp0exp = temp0;
		double tempexp = temp;

		if (functionType == "sig") {
			tempFexp = exp(tempF) / (1 + exp(tempF));
			temp0exp = exp(temp0) / (1 + exp(temp0));
			tempexp = exp(temp) / (1 + exp(temp));
		}
		if (functionType == "log") {
			tempFexp = log(tempF);
			temp0exp = log(temp0);
			tempexp = log(temp);
		}

		if (functionType == "lin") {
			tempFexp = tempF;
			temp0exp = temp0;
			tempexp = temp;
		}

		int currentItMax = static_cast<int>(round((itMax - itMin) / (tempFexp - temp0exp)) * (tempexp - temp0exp) + itMin);
		double currentDelta = (deltaMax - deltaMin) / (tempFexp - temp0exp) * (tempexp - temp0exp) + deltaMin;

		if (functionType == "fix") {
			currentItMax = itMin;
			currentDelta = deltaMin;
		}

		bool melhora = false;
		for (int i = 0; i < currentItMax; i++)
		{
			auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
			if (toc >= FJSSFLE.timeLimit)
				break;
			if (tuning && (toc - ultima_melhora) > 5.0)
				break;

			Gline = G;

			Gline.Perturbate(FJSSFLE);

			double deltaC = ((double)Gline.makespan - (double)G.makespan) / initialMakespan;

			double r = (double)rand() / RAND_MAX;

			if (exp(-deltaC / temp) > r)
			{
				G = Gline;

				if (Gline.makespan < Gbest.makespan) {
#ifdef _DEBUG
					cout << Gbest.makespan << endl;
#endif // _DEBUG

					melhora = true;
					Gbest = Gline;
					auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
					Gbest.timeToGetSolution = toc;
					ultima_melhora = toc;
					ofstream outfile(saidaDetalhada, fstream::app);
					outfile << "SA," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << toc << "," << FJSSFLE.seed << endl;
					outfile.close();
				}
			}
		}

		if (localSearch != "None") {
			G = LocalSearch(FJSSFLE, G, localSearchStrategy, tolerancy, tic);

			if (G.makespan < Gbest.makespan) {
				Gbest = G;
				auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
				Gbest.timeToGetSolution = toc;
				ultima_melhora = toc;
			}
		}
		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
		if (toc >= FJSSFLE.timeLimit)
			break;
		if (tuning && (toc - ultima_melhora) > 5.0)
			break;
		if (melhora)
			itSemMelhora = 0;
		else
			itSemMelhora++;

		temp = max(tempF, temp * currentDelta);
	}

	{
		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "SA," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << toc << "," << FJSSFLE.seed << endl;
		outfile.close();
	}

	return Gbest;
}

SolutionGraph TS(Cenario& FJSSFLE, int listSizeMax, bool tuning, bool onlyCriticalOperations)
{
	auto tic = std::chrono::high_resolution_clock::now();

	SolutionGraph G, Gbest, GSPT, GECT;

	list<Moviment> TabuList;

	GSPT = SPT(FJSSFLE);
	GECT = ECT(FJSSFLE);

	G = (GSPT.makespan < GECT.makespan) ? GSPT : GECT;
	Gbest = G;

	int newListSizeMax = (int)ceil((G.nNodes - 2 + G.nMachines) * listSizeMax / 10);

	string
		saidaDetalhada = "detailed" + FJSSFLE.outputName;

	int nNodes = G.nNodes;
	double ultima_melhora = 0.0;

	{
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "TS," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << 0.0 << "," << FJSSFLE.seed << endl;
		outfile.close();
	}

	while (true) {
		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
		if (toc >= FJSSFLE.timeLimit)
			break;

		if (tuning && (toc - ultima_melhora) > 5.0)
			break;

		for (int i = 0; i < nNodes; i++) {
			isCriticalOpG[i] = G.isACriticalOp[i];
			GTopologicalSorting[i] = G.uu[i];
		}
		int makespanLine = FJSSFLE.BIG_M;
		Moviment movimentoLine;

		for (int i = 0; i < nNodes; i++)
		{
			int v = GTopologicalSorting[i];
			if ((onlyCriticalOperations && !isCriticalOpG[v]) || v == 0 || v == nNodes - 1)
				continue;

			int ellAntigo = G.positions[v],
				kAntigo = G.ff[v],
				MakespanG = G.makespan;

			G.RemoveOperation(FJSSFLE, v);

			for (auto& k : FJSSFLE.F[v - 1])
			{
				int
					gamma_lb = 1,
					gamma_ub = (int)G.Machines[k].size() + 1;

				for (int ell = 0; ell < (int)G.Machines[k].size(); ell++)
				{
					int operation = G.Machines[k][ell];
					if (G.R_in[operation])
						gamma_lb = G.positions[operation] + 1;
					if (G.R_out[operation] && gamma_ub == (int)G.Machines[k].size() + 1)
						gamma_ub = G.positions[operation];
				}

				if (G.makespan >= MakespanG)
					gamma_ub = min(gamma_ub, G.lastCriticalPosition[k]);

				for (int ell = gamma_lb; ell <= gamma_ub; ell++)
				{
					G.InsertOperation(FJSSFLE, v, ell, k);
					Moviment movimento = Moviment(v, k, ell);

					if ((G.makespan < makespanLine && (find(TabuList.begin(), TabuList.end(), movimento) == TabuList.end())) || (G.makespan < Gbest.makespan)) {
						makespanLine = G.makespan;
						movimentoLine = movimento;
					}
					G.RemoveOperation(FJSSFLE, v);
				}
			}
			G.InsertOperation(FJSSFLE, v, ellAntigo, kAntigo);
		}

		auto it = find(TabuList.begin(), TabuList.end(), movimentoLine);
		if (it != TabuList.end())
		{
			TabuList.erase(it);
		}
		else {
			TabuList.push_back(movimentoLine);
		}

		if ((int)TabuList.size() > newListSizeMax)
			TabuList.pop_front();
		
		if (movimentoLine.i == -1 && movimentoLine.ell == -1 && movimentoLine.k == -1) {
			continue;
		}
		else {
			G.RemoveOperation(FJSSFLE, movimentoLine.i);
			G.InsertOperation(FJSSFLE, movimentoLine.i, movimentoLine.ell, movimentoLine.k);
		}

		//atualiza Gbest
		if (makespanLine < Gbest.makespan) {
			Gbest = G;
			auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
			Gbest.timeToGetSolution = toc;
			ultima_melhora = toc;
			ofstream outfile(saidaDetalhada, fstream::app);
			outfile << "TS," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << toc << "," << FJSSFLE.seed << endl;
			outfile.close();

#ifdef _DEBUG
			cout << Gbest.makespan << endl;
#endif
		}
	}

	auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
	if (IMPRIMIR_ARQUIVO) {
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "TS," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << toc << "," << FJSSFLE.seed << endl;
		outfile.close();
	}

	return Gbest;
}

SolutionGraph GRASP(Cenario& FJSSFLE, string localSearchStrategy, double alphaGrasp, double tolerancy, bool tuning, bool onlyCriticalOperations)
{
	auto tic = std::chrono::high_resolution_clock::now();

	SolutionGraph G1, G2, G, Gbest, Gbest_;
	G1 = RandomizedECT(FJSSFLE, 0.0);
	G = ECT(FJSSFLE);
	G2 = RandomizedSPT(FJSSFLE, 0.0);
	G = SPT(FJSSFLE);
	Gbest = G1.makespan < G2.makespan ? G1 : G2;
	double C_max_Best = Gbest.makespan;
	string
		saidaDetalhada = "detailed" + FJSSFLE.outputName;
	size_t pos = saidaDetalhada.find(".txt");
	if (pos != string::npos)
		saidaDetalhada.erase(pos, 4);
	pos = saidaDetalhada.find("instances\\");
	if (pos != string::npos)
		saidaDetalhada.erase(pos, 10);

	int itSemMelhora = 0, iteration = 0;

	{
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "GRASP," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << 0 << "," << FJSSFLE.seed << endl;
		outfile.close();
	}

	int ell = 1;
	double ultima_melhora = 0.0;

	while (true) {
		bool melhora = false;
		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
		if (toc >= FJSSFLE.timeLimit)
			break;
		if (tuning && (toc - ultima_melhora) > 5.0)
			break;

		G1 = RandomizedECT(FJSSFLE, alphaGrasp);
		G2 = RandomizedSPT(FJSSFLE, alphaGrasp);
		G = G1.makespan < G2.makespan ? G1 : G2;

		G = LocalSearch(FJSSFLE, G, localSearchStrategy, tolerancy, tic, onlyCriticalOperations);

		if (G.makespan < C_max_Best) {
			melhora = true;
			Gbest = G;
			C_max_Best = Gbest.makespan;
			auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
			Gbest.timeToGetSolution = toc;
			ultima_melhora = toc;
			ofstream outfile(saidaDetalhada, fstream::app);
			outfile << "GRASP," << FJSSFLE.learningRate << "," << FJSSFLE.inputName << "," << Gbest.makespan << "," << toc << "," << FJSSFLE.seed << endl;
			outfile.close();
			iteration++;
		}
		if (melhora)
			itSemMelhora = 0;
		else
			itSemMelhora++;
		ell++;
	}
	/*if (onlyCriticalOperations)
		Gbest = LocalSearch(FJSSFLE, Gbest, localSearch, localSearchStrategy, 0.0, tic);*/

	return Gbest;
}

SolutionGraph RandomizedECT(Cenario& FJSSFLE, double alphaGrasp) {
	SolutionGraph G = SolutionGraph(FJSSFLE);

	vector<int>
		r_op(G.nNodes, FJSSFLE.BIG_M),
		r_mac(FJSSFLE.nMachines, 0),
		c(G.nNodes, FJSSFLE.BIG_M);
	r_op[0] = 0;
	c[0] = 0;

	vector< int>
		g(FJSSFLE.nMachines, 1);

	vector<bool> isScheduled(G.nNodes, false);
	isScheduled[0] = true;

	set<int> notScheduled;
	for (int i = 1; i < G.nNodes - 1; i++)
		notScheduled.insert(i);

	while (!notScheduled.empty())
	{
		int v_hat = -1, k_hat = -1;
		double
			minCompletionTime = FJSSFLE.BIG_M,
			maxCompletionTime = 0;

		vector<vector<double>> completionTimes(FJSSFLE.nOperations);

		for (int i = 0; i < FJSSFLE.nOperations; i++)
			completionTimes[i] = vector<double>(FJSSFLE.nMachines, FJSSFLE.BIG_M);

		vector< int> ReadyOperations;
		for (auto& v : notScheduled)
		{
			bool isReady = true;
			int readyTime = 0;
			for (auto& i : FJSSFLE.AdjBack[v - 1]) {
				if (!isScheduled[i + 1]) {
					isReady = false;
					break;
				}
				if (c[i + 1] > readyTime)
					readyTime = c[i + 1];
			}
			if (isReady) {
				ReadyOperations.push_back(v);
				r_op[v] = readyTime;
			}
		}
		for (auto& v : notScheduled)
		{
			for (auto& k : FJSSFLE.F[v - 1])
			{
				double realProcessingTime = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[v - 1][k], g[k], FJSSFLE.learningRate);
				double startingTime = max(r_op[v], r_mac[k]);

				completionTimes[v - 1][k] = startingTime + realProcessingTime;
				if (completionTimes[v - 1][k] < minCompletionTime)
					minCompletionTime = completionTimes[v - 1][k];

				if (completionTimes[v - 1][k] > maxCompletionTime)
					maxCompletionTime = completionTimes[v - 1][k];
			}
		}

		vector<Moviment> RCL;
		for (auto& v : ReadyOperations)
		{
			for (auto& k : FJSSFLE.F[v - 1]) {
				if (completionTimes[v - 1][k] <= minCompletionTime + alphaGrasp * (maxCompletionTime - minCompletionTime))
					RCL.push_back(Moviment(v, k));
			}
		}
		int chosenMoviment = alphaGrasp == 0.0 ? 0 : rand() % ((int)RCL.size());
		v_hat = RCL[chosenMoviment].i;
		k_hat = RCL[chosenMoviment].k;

		G.weight[v_hat] = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[v_hat - 1][k_hat], g[k_hat], FJSSFLE.learningRate);
		G.ff[v_hat] = k_hat;
		c[v_hat] = max(r_op[v_hat], r_mac[k_hat]) + G.weight[v_hat];
		r_mac[k_hat] = c[v_hat];
		G.positions[v_hat] = g[k_hat];
		G.st[v_hat] = max(r_op[v_hat], r_mac[k_hat]);
		g[k_hat]++;

		if (!G.Machines[k_hat].empty())
		{
			G.AdjList[G.Machines[k_hat].back()].push_back(SolutionNode(v_hat, k_hat, G.weight[v_hat]));
		}
		G.Machines[k_hat].push_back(v_hat);
		notScheduled.erase(v_hat);
		isScheduled[v_hat] = true;
		//cout << "(" << v_hat << "," << k_hat << ")" << endl;
	}
	G.CriticalPath(FJSSFLE);
	return G;
}

SolutionGraph RandomizedSPT(Cenario& FJSSFLE, double alphaGrasp) {
	SolutionGraph G = SolutionGraph(FJSSFLE);

	vector<int>
		r_op(G.nNodes, FJSSFLE.BIG_M),
		r_mac(FJSSFLE.nMachines, 0),
		c(G.nNodes, FJSSFLE.BIG_M);
	r_op[0] = 0;
	c[0] = 0;

	vector< int>
		g(FJSSFLE.nMachines, 1);

	vector<bool> isScheduled(G.nNodes, false);
	isScheduled[0] = true;

	set<int> notScheduled;
	for (int i = 1; i < G.nNodes - 1; i++)
		notScheduled.insert(i);

	while (!notScheduled.empty())
	{
		int v_hat = -1, k_hat = -1;
		int
			minProcessingTime = FJSSFLE.BIG_M,
			maxProcessingTime = 0;

		vector<int> ReadyOperations;

		int r_min = FJSSFLE.BIG_M;
		for (auto& v : notScheduled)
		{
			bool isReady = true;
			int readyTime = 0;
			for (auto& i : FJSSFLE.AdjBack[v - 1]) {
				if (!isScheduled[i + 1]) {
					isReady = false;
					break;
				}
				if (c[i + 1] > readyTime)
					readyTime = c[i + 1];
			}
			if (isReady) {
				r_op[v] = readyTime;
				ReadyOperations.push_back(v);

				for (auto& k : FJSSFLE.F[v - 1])
				{
					if (r_min > max(r_op[v], r_mac[k]))
						r_min = max(r_op[v], r_mac[k]);
				}
			}
		}

		vector<Moviment> E;
		for (auto& v : ReadyOperations)
		{
			for (auto& k : FJSSFLE.F[v - 1])
			{
				if (r_min == max(r_op[v], r_mac[k]))
				{
					E.push_back(Moviment(v, k));
					int realProcessingTime = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[v - 1][k], g[k], FJSSFLE.learningRate);

					if (realProcessingTime < minProcessingTime)
					{
						minProcessingTime = realProcessingTime;
					}
					if (realProcessingTime > maxProcessingTime)
					{
						maxProcessingTime = realProcessingTime;
					}
				}
			}
		}

		vector<Moviment> RCL;
		for (auto& moviment : E)
		{
			int
				v = moviment.i,
				k = moviment.k;

			double realProcessingTime = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[v - 1][k], g[k], FJSSFLE.learningRate);

			if (realProcessingTime <= minProcessingTime + alphaGrasp * (maxProcessingTime - minProcessingTime))
				RCL.push_back(Moviment(v, k));
		}
		int chosenMoviment = alphaGrasp == 0.0 ? 0 : rand() % ((int)RCL.size());
		v_hat = RCL[chosenMoviment].i;
		k_hat = RCL[chosenMoviment].k;

		G.weight[v_hat] = FJSSFLE.learningFunction(FJSSFLE.ProcessingTimes[v_hat - 1][k_hat], g[k_hat], FJSSFLE.learningRate);
		G.ff[v_hat] = k_hat;
		c[v_hat] = r_min + G.weight[v_hat];
		r_mac[k_hat] = c[v_hat];
		G.positions[v_hat] = g[k_hat];
		G.st[v_hat] = r_min;
		g[k_hat]++;

		if (!G.Machines[k_hat].empty())
		{
			G.AdjList[G.Machines[k_hat].back()].push_back(SolutionNode(v_hat, k_hat, G.weight[v_hat]));
		}
		G.Machines[k_hat].push_back(v_hat);
		notScheduled.erase(v_hat);
		isScheduled[v_hat] = true;
		//cout << "(" << v_hat << "," << k_hat << ")" << endl;
	}
	G.CriticalPath(FJSSFLE);
	return G;
}