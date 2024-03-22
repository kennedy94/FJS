#include "Tayebi.h"

void evaluateChromosome(FJS& P, Chromosome& S) {
	P.nAvaliacoesMakespan++;
	int fitness = 0;

	vector<int>
		time_machine(P.nMac, 0),
		ct(P.nOp, 0);
	vector<int>
		positions(P.nMac, 0);

	for (int i = 0; i < P.nOp; i++) {
		int maximo_ct = 0;

		for (auto& prec : P.AdjBack[S.operations[i]]) {
			if (ct[prec] > maximo_ct)
				maximo_ct = ct[prec];
		}
		positions[S.machines[i]]++;
		int actual_processing_time =
			learningFunction(P.processingTimesMatrix[S.operations[i]][S.machines[i]], positions[S.machines[i]], P.alpha);

		time_machine[S.machines[i]] = max(maximo_ct, time_machine[S.machines[i]]) + actual_processing_time;
		ct[S.operations[i]] = time_machine[S.machines[i]];
		fitness = max(fitness, ct[S.operations[i]]);
	}

	S.fitness = fitness;
}

void repairChromosome(FJS& P, Chromosome& S) {
	vector<vector< int>>
		JobsIdx(P.nJobs);

	vector< int> ff(P.nOp);

	for (int i = 0; i < P.nOp; i++) {
		int operacao = S.operations[i];
		ff[operacao] = S.machines[i];
		JobsIdx[P.operationsVector[operacao].job].push_back(i);
	}
	//O(n) amortizado

	for (int j = 0; j < P.nJobs; j++) {
		for (int i = 0; i < (int)JobsIdx[j].size(); i++) {
			S.operations[JobsIdx[j][i]] = P.Jobs[j][i];
			S.machines[JobsIdx[j][i]] = ff[P.Jobs[j][i]];
		}
	}
}

vector<Chromosome> generatePopulation(FJS& P, int popsize)
{
	mt19937 generator((int)time(NULL));

	vector< int>
		permutation(P.nOp),
		machines(P.nOp, 0);

	vector<Chromosome> Population(popsize);

	for (int i = 0; i < popsize; i++) {
		for (int j = 0; j < (int)permutation.size(); ++j)
			permutation[j] = j;

		random_shuffle(permutation.begin(), permutation.end());

		Chromosome S;
		S.operations = permutation;
		for (int j = 0; j < (int)permutation.size(); ++j) {
			int op = S.operations[j];
			uniform_int_distribution<int> distribution(0, P.operationsVector[op].Fsize - 1);
			int number = distribution(generator);
			machines[j] = P.operationsVector[op].F[number];
		}

		S.machines = machines;

		repairChromosome(P, S);
		evaluateChromosome(P, S);
		Population[i] = S;
	}

	return Population;
}

Chromosome neighbourhood(FJS& P, Chromosome& S, int movimento)
{
	int i = rand() % ((int)S.machines.size()),
		j;
	do
	{
		j = rand() % ((int)S.machines.size());
	} while (i == j);

	if (i > j) {
		int aux = i;
		i = j;
		j = aux;
	}

	vector<string> scenarios = { "One", "Two", "Three" };
	int moeda = rand() % 3;

	switch (movimento)
	{
	case 1:
		return swap(P, S, scenarios[moeda], i, j);
	case 2:
		return inversion(P, S, scenarios[moeda], i, j);
	case 3:
		return insertion(P, S, scenarios[moeda], i, j);
	case 4:
		return threepoint(P, S, scenarios[moeda]);
	case 5:
		return fourpoint(P, S, scenarios[moeda]);
	default:
		return S;
	}
}

Chromosome insertion(FJS& P, Chromosome& S, string scenario, int i, int j)
{
	Chromosome vizinho = S;
	if (scenario == "One" || scenario == "Three") {
		vizinho.operations[i + 1] = S.operations[j];

		for (int it = i + 2; it <= j; it++) {
			vizinho.operations[it] = S.operations[it - 1];
			if (it == (int)S.machines.size() - 1)
				break;
		}
	}
	if (scenario == "Two" || scenario == "Three") {
		vizinho.machines[i + 1] = S.machines[j];
		if (find(P.operationsVector[vizinho.operations[i + 1]].F.begin(), P.operationsVector[vizinho.operations[i + 1]].F.end(), vizinho.machines[i + 1])
			== P.operationsVector[vizinho.operations[i + 1]].F.end()) {
			int random_machine = rand() % ((int)P.operationsVector[vizinho.operations[i + 1]].F.size());
			vizinho.machines[i + 1] = P.operationsVector[vizinho.operations[i + 1]].F[random_machine];
		}

		for (int it = i + 2; it <= j; it++) {
			vizinho.machines[it] = S.machines[it - 1];

			if (find(P.operationsVector[vizinho.operations[it]].F.begin(), P.operationsVector[vizinho.operations[it]].F.end(), vizinho.machines[it])
				== P.operationsVector[vizinho.operations[it]].F.end()) {
				int random_machine = rand() % ((int)P.operationsVector[vizinho.operations[it]].F.size());
				vizinho.machines[it] = P.operationsVector[vizinho.operations[it]].F[random_machine];
			}

			if (it == (int)S.machines.size() - 1)
				break;
		}
	}

	return vizinho;
}

Chromosome swap(FJS& P, Chromosome& S, string scenario, int i, int j)
{
	Chromosome vizinho = S;
	if (scenario == "One" || scenario == "Three") {
		vizinho.operations[i] = S.operations[j];
		vizinho.operations[j] = S.operations[i];
	}
	if (scenario == "Two" || scenario == "Three") {
		vizinho.machines[i] = S.machines[j];
		vizinho.machines[j] = S.machines[i];

		if (find(P.operationsVector[vizinho.operations[i]].F.begin(), P.operationsVector[vizinho.operations[i]].F.end(), vizinho.machines[i])
			== P.operationsVector[vizinho.operations[i]].F.end()) {
			int random_machine = rand() % ((int)P.operationsVector[vizinho.operations[i]].F.size());
			vizinho.machines[i] = P.operationsVector[vizinho.operations[i]].F[random_machine];
		}

		if (find(P.operationsVector[vizinho.operations[j]].F.begin(), P.operationsVector[vizinho.operations[j]].F.end(), vizinho.machines[j])
			== P.operationsVector[vizinho.operations[j]].F.end()) {
			int random_machine = rand() % ((int)P.operationsVector[vizinho.operations[j]].F.size());
			vizinho.machines[j] = P.operationsVector[vizinho.operations[j]].F[random_machine];
		}
	}

	return vizinho;
}

void crossover(FJS& P, vector<Chromosome>& pop2, Chromosome& P1, Chromosome& P2, string mode, string scenario) {
	if (P1.operations.size() != P2.operations.size() || P1.machines.size() != P2.machines.size()) {
		cerr << "Erro no crossover, dimensões dos pais diferentes" << endl;
		exit(-1);
	}
	int n = (int)P1.operations.size();
	Chromosome
		child1 = P1,
		child2 = P2;

	if (mode == "Uniform") {
		for (int i = 0; i < n; i++) {
			double moeda = (double)rand() / RAND_MAX;
			vector<int> mapping(P.nOp, -1), mapping2(P.nOp, -1);

			if (moeda < 0.5) {
				if (scenario == "One" || scenario == "Three") {
					mapping[P1.operations[i]] = P2.operations[i];
					mapping2[P2.operations[i]] = P1.operations[i];

					child1.operations[i] = P2.operations[i];
					child2.operations[i] = P1.operations[i];
				}
				if (scenario == "Two" || scenario == "Three") {
					//se a máquina não processa a operação
					if (find(P.operationsVector[child1.operations[i]].F.begin(), P.operationsVector[child1.operations[i]].F.end(), P2.machines[i]) == P.operationsVector[child1.operations[i]].F.end())
					{
						int random_machine;
						do
						{
							random_machine = rand() % P.operationsVector[child1.operations[i]].F.size();
							random_machine = P.operationsVector[child1.operations[i]].F[random_machine];
						} while (random_machine == child1.machines[i] && P.operationsVector[child1.operations[i]].F.size() > 1);
						child1.machines[i] = random_machine;
					}
					else
						child1.machines[i] = P2.machines[i];

					if (find(P.operationsVector[child2.operations[i]].F.begin(), P.operationsVector[child2.operations[i]].F.end(), P1.machines[i]) == P.operationsVector[child2.operations[i]].F.end()) {
						int random_machine;
						do
						{
							random_machine = rand() % P.operationsVector[child2.operations[i]].F.size();
							random_machine = P.operationsVector[child2.operations[i]].F[random_machine];
						} while (random_machine == child2.machines[i] && P.operationsVector[child2.operations[i]].F.size() > 1);
						child2.machines[i] = random_machine;
					}
					else
						child2.machines[i] = P1.machines[i];
				}
			}
			vector< int> frequency(P.nOp, 0), frequency2(P.nOp, 0);
			for (int i = 0; i < P.nOp; i++) {
				frequency[child1.operations[i]]++;
				frequency2[child2.operations[i]]++;
			}

			for (int i = 0; i < P.nOp; i++) {
				if (frequency[child1.operations[i]] == 2) {
					frequency[child1.operations[i]]--;
					int aux = child1.operations[i];
					while (aux != -1) {
						child1.operations[i] = mapping2[child1.operations[i]];
						aux = mapping2[child1.operations[i]];
					}
				}

				if (frequency2[child2.operations[i]] == 2) {
					frequency2[child2.operations[i]]--;
					int aux = child2.operations[i];
					while (aux != -1) {
						child2.operations[i] = mapping[child2.operations[i]];
						aux = mapping[child2.operations[i]];
					}
				}
			}
		}
	}
	if (mode == "OnePoint") {
		int
			OperationsPoint = rand() % (n - 1),
			MachinesPoint = rand() % (n - 1);

		vector<int> mapping(P.nOp, -1), mapping2(P.nOp, -1);

		for (int i = 0; i < n; i++) {
			if (scenario == "One" || scenario == "Three") {
				if (i > OperationsPoint) {
					mapping[P1.operations[i]] = P2.operations[i];
					mapping2[P2.operations[i]] = P1.operations[i];

					child1.operations[i] = P2.operations[i];
					child2.operations[i] = P1.operations[i];
				}
			}
			if (scenario == "Two" || scenario == "Three") {
				if (i > MachinesPoint) {
					//se a máquina não processa a operação
					if (find(P.operationsVector[child1.operations[i]].F.begin(), P.operationsVector[child1.operations[i]].F.end(), P2.machines[i]) == P.operationsVector[child1.operations[i]].F.end())
					{
						int random_machine;
						do
						{
							random_machine = rand() % P.operationsVector[child1.operations[i]].F.size();
							random_machine = P.operationsVector[child1.operations[i]].F[random_machine];
						} while (random_machine == child1.machines[i] && P.operationsVector[child1.operations[i]].F.size() > 1);
						child1.machines[i] = random_machine;
					}
					else
						child1.machines[i] = P2.machines[i];

					if (find(P.operationsVector[child2.operations[i]].F.begin(), P.operationsVector[child2.operations[i]].F.end(), P1.machines[i]) == P.operationsVector[child2.operations[i]].F.end()) {
						int random_machine;
						do
						{
							random_machine = rand() % P.operationsVector[child2.operations[i]].F.size();
							random_machine = P.operationsVector[child2.operations[i]].F[random_machine];
						} while (random_machine == child2.machines[i] && P.operationsVector[child2.operations[i]].F.size() > 1);
						child2.machines[i] = random_machine;
					}
					else
						child2.machines[i] = P1.machines[i];
				}
			}
		}
		vector< int> frequency(P.nOp, 0), frequency2(P.nOp, 0);
		for (int i = 0; i < P.nOp; i++) {
			frequency[child1.operations[i]]++;
			frequency2[child2.operations[i]]++;
		}

		for (int i = 0; i < P.nOp; i++) {
			if (frequency[child1.operations[i]] == 2) {
				frequency[child1.operations[i]]--;
				int aux = child1.operations[i];
				while (aux != -1) {
					child1.operations[i] = mapping2[child1.operations[i]];
					aux = mapping2[child1.operations[i]];
				}
			}

			if (frequency2[child2.operations[i]] == 2) {
				frequency2[child2.operations[i]]--;
				int aux = child2.operations[i];
				while (aux != -1) {
					child2.operations[i] = mapping[child2.operations[i]];
					aux = mapping[child2.operations[i]];
				}
			}
		}
	}

	repairChromosome(P, child1);
	evaluateChromosome(P, child1);
	pop2.push_back(child1);

	repairChromosome(P, child2);
	evaluateChromosome(P, child2);
	pop2.push_back(child2);
}

bool isFeasible(FJS& P, Chromosome S)
{
	vector< int> countJob(P.nJobs, 0);

	for (int i = 0; i < P.nOp; i++) {
		int operacao = S.operations[i];
		if (operacao != P.Jobs[P.operationsVector[operacao].job][countJob[P.operationsVector[operacao].job]])
			return false;

		if (find(P.operationsVector[operacao].F.begin(), P.operationsVector[operacao].F.end(), S.machines[i]) == P.operationsVector[operacao].F.end())
			return false;

		countJob[P.operationsVector[operacao].job]++;
	}

	return true;
}

Chromosome inversion(FJS& P, Chromosome& S, string scenario, int i, int j)
{
	Chromosome vizinho = S;
	if (scenario == "One" || scenario == "Three") {
		for (int it = 0; it <= j - i; it++) {
			vizinho.operations[i + it] = S.operations[j - it];
		}
	}
	if (scenario == "Two" || scenario == "Three") {
		for (int it = 0; it <= j - i; it++) {
			vizinho.machines[i + it] = S.machines[j - it];

			if (find(P.operationsVector[vizinho.operations[i + it]].F.begin(), P.operationsVector[vizinho.operations[i + it]].F.end(), vizinho.machines[i + it])
				== P.operationsVector[vizinho.operations[i + it]].F.end()) {
				int random_machine = rand() % P.operationsVector[vizinho.operations[i + it]].F.size();
				vizinho.machines[i + it] = P.operationsVector[vizinho.operations[i + it]].F[random_machine];
			}
		}
	}

	return vizinho;
}

Chromosome threepoint(FJS& P, Chromosome& S, string scenario)
{
	vector< int> idx(3, 0);
	int n = (int)S.machines.size();

	idx[0] = rand() % n;
	do
	{
		idx[1] = rand() % n;
		do
		{
			idx[2] = rand() % n;
		} while (idx[1] == idx[2] || idx[0] == idx[2]);
	} while (idx[0] == idx[1]);

	sort(idx.begin(), idx.end());
	vector< int> idxO = idx;

	Chromosome
		vizinho = S;

	while (next_permutation(idx.begin(), idx.end())) {
		Chromosome aux = S;
		if (scenario == "One" || scenario == "Three") {
			for (int i = 0; i < 3; i++)
				aux.operations[idx[i]] = S.operations[idxO[i]];
		}
		if (scenario == "Two" || scenario == "Three") {
			for (int i = 0; i < 3; i++) {
				aux.machines[idx[i]] = S.machines[idxO[i]];
				if (find(P.operationsVector[vizinho.operations[idx[i]]].F.begin(), P.operationsVector[vizinho.operations[idx[i]]].F.end(), vizinho.machines[idx[i]])
					== P.operationsVector[vizinho.operations[idx[i]]].F.end()) {
					int random_machine = rand() % P.operationsVector[vizinho.operations[idx[i]]].F.size();
					vizinho.machines[idx[i]] = P.operationsVector[vizinho.operations[idx[i]]].F[random_machine];
				}
			}
		}
		repairChromosome(P, aux);
		evaluateChromosome(P, aux);
		if (aux.fitness < vizinho.fitness)
			vizinho = aux;
	}

	return vizinho;
}

Chromosome fourpoint(FJS& P, Chromosome& S, string scenario)
{
	vector< int> idx(4, 0);
	int n = (int)S.machines.size();

	idx[0] = rand() % n;
	do
	{
		idx[1] = rand() % n;
		do
		{
			idx[2] = rand() % n;
			do
			{
				idx[3] = rand() % n;
			} while (idx[2] == idx[3] || idx[0] == idx[3] || idx[1] == idx[3]);
		} while (idx[1] == idx[2] || idx[0] == idx[2]);
	} while (idx[0] == idx[1]);

	sort(idx.begin(), idx.end());
	vector< int> idxO = idx;

	Chromosome
		vizinho = S;

	while (next_permutation(idx.begin(), idx.end())) {
		Chromosome aux = S;
		if (scenario == "One" || scenario == "Three") {
			for (int i = 0; i < 4; i++)
				aux.operations[idx[i]] = S.operations[idxO[i]];
		}
		if (scenario == "Two" || scenario == "Three") {
			for (int i = 0; i < 4; i++) {
				aux.machines[idx[i]] = S.machines[idxO[i]];
				if (find(P.operationsVector[vizinho.operations[idx[i]]].F.begin(), P.operationsVector[vizinho.operations[idx[i]]].F.end(), vizinho.machines[idx[i]])
					== P.operationsVector[vizinho.operations[idx[i]]].F.end()) {
					int random_machine = rand() % P.operationsVector[vizinho.operations[idx[i]]].F.size();
					vizinho.machines[idx[i]] = P.operationsVector[vizinho.operations[idx[i]]].F[random_machine];
				}
			}
		}
		repairChromosome(P, aux);
		evaluateChromosome(P, aux);
		if (aux.fitness < vizinho.fitness)
			vizinho = aux;
	}

	return vizinho;
}

vector<Chromosome> xSelections(vector<Chromosome>& popu, int M)
{
	vector< int> selected;
	vector<Chromosome> selectedChr;

	for (int i = 0; i < min(M, (int)popu.size()); i++) {
		int rand_idx;
		do
			rand_idx = rand() % popu.size();
		while (find(selected.begin(), selected.end(), rand_idx) != selected.end());
		selected.push_back(rand_idx);
		selectedChr.push_back(popu[rand_idx]);
	}

	sort(selectedChr.begin(), selectedChr.end());

	return selectedChr;
}

Chromosome TayebiGA(FJS& P, int npop, int M,
	double Pc, double Pns, double Paf, int Lmax)
{
	vector<Chromosome> pop1 = generatePopulation(P, npop);
	sort(pop1.begin(), pop1.end());

	auto t1 = std::chrono::high_resolution_clock::now();
	auto toc = chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count() / 1000.0;

	string
		saidaDetalhada = "detailed" + P.outputName;
	size_t pos = saidaDetalhada.find(".txt");
	if (pos != string::npos)
		saidaDetalhada.erase(pos, 4);
	pos = saidaDetalhada.find("instances\\");
	if (pos != string::npos)
		saidaDetalhada.erase(pos, 10);

	int iter_sem_melhora = 0;
	double best_makespan_antigo = 0.0;

	{
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "Tayebi," << P.alpha << "," << P.instanceName << "," << pop1.front().fitness << "," << 0.0 << "," << P.SEED << endl;
		outfile.close();
	}

	//for (int i = 0; i < MaxIt; i++) {
	while (true) {
		vector<Chromosome> pop2;
		for (int j = 0; j < round((double)Pc * npop / 2.0); j++)
		{
			vector<Chromosome> Selected = xSelections(pop1, M);

			string scenario,
				mode;

			double moeda = (double)rand() / RAND_MAX;
			if (moeda < 0.33333)
				scenario = "One";
			else
			{
				if (moeda < 0.66666)
					scenario = "Two";
				else
					scenario = "Three";
			}
			moeda = (double)rand() / RAND_MAX;
			if (moeda < 0.5)
				mode = "Uniform";
			else
				mode = "OnePoint";
			crossover(P, pop2, Selected[0], Selected[1], mode, scenario);
		}

		vector<Chromosome> pop3;
		for (int j = 0; j < round((double)Pns * npop / 2.0); j++) {
			vector<Chromosome> Selected = xSelections(pop2, M);
			Chromosome BEST = Selected.front(), X = Selected.front();

			for (int L = 0; L < Lmax; L++) {
				for (int K = 1; K <= 5; K++) {
					Chromosome LK = neighbourhood(P, X, K);
					repairChromosome(P, LK);
					evaluateChromosome(P, LK);
					if (LK.fitness < BEST.fitness) {
						BEST = LK;
						pop3.push_back(BEST);
						break;
					}
				}
			}
		}

		for (auto& child : pop2) {
			pop1.push_back(child);
		}
		for (auto& child : pop3) {
			pop1.push_back(child);
		}

		vector<Chromosome> unique_pop;

		best_makespan_antigo = pop1.front().fitness;

		for (auto& individual : pop1) {
			if (find(unique_pop.begin(), unique_pop.end(), individual) == unique_pop.end()) {
				unique_pop.push_back(individual);
			}
		}
		sort(unique_pop.begin(), unique_pop.end());
		pop1.clear();
		for (int it = 0; it < (int)(Paf * (double)unique_pop.size()); it++) {
			pop1.push_back(unique_pop[it]);
			if ((int)pop1.size() >= npop)
				break;
		}

		int max = (int)unique_pop.size() - 1,
			min = (int)pop1.size();

		while ((int)pop1.size() < npop) {
			int roleta_uniforme = rand() % (max - min + 1) + min;
			pop1.push_back(unique_pop[roleta_uniforme]);
		}

		sort(pop1.begin(), pop1.end());

		toc = chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count() / 1000.0;
		if (best_makespan_antigo > pop1.front().fitness)
		{
			ofstream outfile(saidaDetalhada, fstream::app);
			outfile << "Tayebi," << P.alpha << "," << P.instanceName << "," << pop1.front().fitness << "," << toc << "," << P.SEED << endl;
			outfile.close();

			iter_sem_melhora = 0;
			pop1.front().time_to_find = toc;
		}
		else {
			iter_sem_melhora++;
		}

		if (toc > P.TimeLimit)
			break;
	}

	{
		ofstream outfile(saidaDetalhada, fstream::app);
		outfile << "Tayebi," << P.alpha << "," << P.instanceName << "," << pop1.front().fitness << "," << toc << "," << P.SEED << endl;
		outfile.close();
	}

	if (isFeasible(P, pop1.front()))
		cout << "Viavel, tudo okay!" << endl;
	else
		cout << "Inviavel :(" << endl;

	return pop1.front();
}