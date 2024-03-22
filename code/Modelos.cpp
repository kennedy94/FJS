#include "Modelos.h"

bool ModeloFeasibilityTest(const SolutionGraph& G, Cenario& Problem) {
	vector<vector<vector<bool>>> x(Problem.nOperations);
	vector<int> s(Problem.nOperations, 0);
	vector<int> p_(Problem.nOperations, 0);
	vector<vector<int>> h(Problem.nMachines);

	for (int i = 0; i < Problem.nOperations; i++)
	{
		s[i] = G.st[i + 1];
		p_[i] = G.weight[i + 1];
	}

	for (IloInt k = 0; k < Problem.nMachines; k++) {
		h[k] = vector<int>((int)Problem.I[k].size(), Problem.BIG_M);
		for (int r = 0; r < (int)G.Machines[k].size(); r++)
		{
			h[k][r] = G.st[G.Machines[k][r]];
		}
	}

	for (int i = 0; i < Problem.nOperations; i++) {
		x[i] = vector<vector<bool>>(Problem.nMachines);
		for (auto& k : Problem.OperationSet[i].machineSet) {
			x[i][k] = vector<bool>((int)Problem.I[k].size(), false);
		}
	}
	for (int i = 0; i < Problem.nOperations; i++)
		x[i][G.ff[i + 1]][G.positions[i + 1] - 1] = true;

	//Restrições
	for (int i = 0; i < Problem.nOperations; i++) {
		int sum = 0;
		for (auto k : Problem.OperationSet[i].machineSet) {
			for (int r = 0; r < (int)Problem.I[k].size(); r++) {
				sum += x[i][k][r];
			}
		}
		if (sum != 1)
			return false;
	}

	for (IloInt i = 0; i < Problem.nOperations; i++) {
		int sum = 0;
		for (auto k : Problem.OperationSet[i].machineSet) {
			for (int r = 0; r < (int)Problem.I[k].size(); r++) {
				sum += x[i][k][r];
			}
		}
		if (sum != 1)
			return false;

		sum = 0;

		for (auto& k : Problem.OperationSet[i].machineSet) {
			for (int r = 0; r < (int)Problem.I[k].size(); r++) {
				sum += x[i][k][r] * Problem.learningFunction(Problem.ProcessingTimes[i][k], r + 1, Problem.learningRate);
			}
		}

		if (p_[i] != sum)
			return false;
		if (s[i] + p_[i] > G.makespan)
			return false;
	}
	return true;
}

void modelo_MILP_posicional_CPLEX(Cenario& Problem) {
	IloEnv env;
	IloModel model(env, "MILP from Section 3.1");
	int nIntegerVariables = 0,
		nContinousVariables = 0,
		nConstraints = 0;

	/* Declaring model variables */
	IloNumVar Cmax(env, 0.0, IloInfinity, ILOFLOAT);
	nContinousVariables++;
	IloNumVarArray s(env, Problem.nOperations, 0.0, IloInfinity, ILOFLOAT);
	nContinousVariables += Problem.nOperations;
	IloNumVarArray p_(env, Problem.nOperations, 0.0, IloInfinity, ILOFLOAT);
	nContinousVariables += Problem.nOperations;

	IloArray<IloArray<IloBoolVarArray>> x(env, Problem.nOperations);
	IloArray<IloNumVarArray> h(env, Problem.nMachines);
	nContinousVariables += Problem.nMachines;

	for (IloInt i = 0; i < Problem.nOperations; i++) {
		x[i] = IloArray<IloBoolVarArray>(env, (IloInt)Problem.nMachines);
		for (auto &k : Problem.OperationSet[i].machineSet) {
			x[i][k] = IloBoolVarArray(env, (IloInt)Problem.I[k].size());
			nIntegerVariables += (int)Problem.I[k].size();
		}
	}

	for (IloInt k = 0; k < Problem.nMachines; k++) {
		h[k] = IloNumVarArray(env, (IloInt)Problem.I[k].size(), 0.0, IloInfinity, ILOFLOAT);
		nContinousVariables += (int)Problem.I[k].size();
	}

	model.add(IloMinimize(env, Cmax));

	for (IloInt i = 0; i < Problem.nOperations; i++) {
		IloExpr sum(env);
		for (auto &k : Problem.OperationSet[i].machineSet) {
			for (IloInt r = 0; r < (IloInt)Problem.I[k].size(); r++) {
				sum += x[i][k][r];
			}
		}
		model.add(sum == 1); //(2)
		nConstraints++;
		sum.clear();

		for (IloInt k = 0; k < (IloInt)Problem.OperationSet[i].nMachines; k++) {
			for (IloInt r = 0; r < (IloInt)Problem.I[Problem.OperationSet[i].machineSet[k]].size(); r++) {
				int m = Problem.OperationSet[i].machineSet[k];
				sum += x[i][m][r] * Problem.learningFunction(Problem.ProcessingTimes[i][m], (int)r + 1, Problem.learningRate);
			}
		}

		model.add(p_[i] == sum).setName("processing time");
		nConstraints++;
		sum.clear(); //(5)

		model.add(s[i] + p_[i] <= Cmax); //(5)
	}

	for (IloInt k = 0; k < (IloInt)Problem.nMachines; k++) {
		for (IloInt r = 0; r < (IloInt)Problem.I[k].size(); r++) {
			IloExpr sum(env), sum1(env);
			for (auto& i : Problem.I[k]) {
				sum += x[i][k][r];
				if (r != ((int)Problem.I[k].size() - 1))
					sum1 += x[i][k][r + 1];
			}
			model.add(sum <= 1); //(3)
			nConstraints++;
			model.add(sum1 <= sum); //(4)
			nConstraints++;
			sum.clear();
			sum1.clear();
		}
	}

	for (IloInt a = 0; a < Problem.nArcs; a++) {
		model.add(s[Problem.PrecedenceArcSet[a].source] + p_[Problem.PrecedenceArcSet[a].source] <= s[Problem.PrecedenceArcSet[a].tail]).setName("Precedence"); //(8)
		nConstraints++;
	}

	for (IloInt k = 0; k < (IloInt)Problem.nMachines; k++) {
		for (IloInt r = 0; r < (IloInt)Problem.I[k].size(); r++) {
			IloExpr sum(env);
			for (auto i : Problem.I[k]) {
				sum += x[i][k][r] * Problem.learningFunction(Problem.ProcessingTimes[i][k], (int)r + 1, Problem.learningRate);;
			}
			if (r == (IloInt)Problem.I[k].size() - 1) {
				model.add(h[k][r] + sum <= Cmax).setName("h"); //(7)
				nConstraints++;
			}
			else {
				model.add(h[k][r] + sum <= h[k][r + 1]).setName("h"); //(6)
				nConstraints++;
			}
			sum.clear();
		}
	}

	for (IloInt i = 0; i < Problem.nOperations; i++) {
		for (IloInt j = 0; j < Problem.nOperations; j++) {
			if (i != j) {
				for (auto k : Problem.OperationSet[i].machineSet) {
					bool found = (Problem.OperationSet[j].machineSet.end() != find(Problem.OperationSet[j].machineSet.begin(), Problem.OperationSet[j].machineSet.end(), k));

					if (found) {
						for (IloInt r = 0; r < (IloInt)Problem.I[k].size() - 1; r++) {
							IloExpr sum(env);
							for (IloInt rr = r + 1; rr < (IloInt)Problem.I[k].size(); rr++)
								sum += x[j][k][rr];

							model.add(s[i] + p_[i] - Problem.BIG_M * (2 - x[i][k][r] - sum) <= s[j]).setName("s"); //(9)
							nConstraints++;

							sum.end();
						}
					}
				}
			}
		}
	}

	for (IloInt i = 0; i < Problem.nOperations; i++)
	{
		for (auto k : Problem.OperationSet[i].machineSet) {
			for (IloInt r = 0; r < (IloInt)Problem.I[k].size(); r++) {
				//model.add(s[i] - M * (1 - x[i][r][k]) <= h[r][k]);

				//model.add(s[i] + M * (1 - x[i][r][k]) >= h[r][k]);

				model.add(h[k][r] - Problem.BIG_M * (1 - x[i][k][r]) <= s[i]).setName("s"); //(11)
				nConstraints++;

				model.add(s[i] <= h[k][r] + Problem.BIG_M * (1 - x[i][k][r])).setName("s"); //(10)
				nConstraints++;
			}
		}
	}

	IloCplex cplex(model);
	ofstream instancesInfo("instancesInfo.csv", fstream::app);
	instancesInfo
		<< "CPLEX,"
		<< Problem.inputName << ","
		<< nIntegerVariables << ","
		<< nContinousVariables << ","
		<< nConstraints << endl;
	instancesInfo.close();
	//return;

	//cplex.setOut("saida.out");
	cplex.setParam(IloCplex::Param::TimeLimit, (IloNum)Problem.timeLimit);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
//#ifdef _DEBUG
//	cplex.setParam(IloCplex::Param::MIP::Display, 1);
//#endif // _DEBUG

	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1 - 1e-06);
	cplex.setParam(IloCplex::Param::WorkMem, 4000);
	cplex.setParam(IloCplex::Param::MIP::Strategy::File, 3);
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	//cplex.setParam(IloCplex::ClockType, 1); /* CPU time */
	//cplex.setParam(IloCplex::WorkMem, 2048); /* Memory limit */

	bool warmStart = true;
	if (warmStart) {
		SolutionGraph G, GSPT, GECT;
		GSPT = SPT(Problem);
		GECT = ECT(Problem);

		G = (GSPT.makespan < GECT.makespan) ? GSPT : GECT;
		//G = LocalSearch(Problem, G, "Best", 0, std::chrono::high_resolution_clock::now(), 0);
		//ModeloFeasibilityTest(G, Problem);
		IloSolution startsol(env);
		{
			IloNumVarArray startVar(env);
			IloNumArray startVal(env);
			for (IloInt i = 0; i < Problem.nOperations; i++) {
				startVar.add(x[i][G.ff[i + 1]][G.positions[i + 1] - 1]);
				startVal.add(1);

				startVar.add(s[i]);
				startVal.add(G.st[i + 1]);
			}
			cplex.addMIPStart(startVar, startVal);
		}
	}

#if _DEBUG
	cplex.exportModel("modelo_alternativo.lp");
#endif

	IloNum elapsed;
	elapsed = cplex.getTime();

	if (!cplex.solve()) {
		elapsed = cplex.getTime() - elapsed;
		ofstream output(Problem.outputName, fstream::app);
		output
			<< "CPLEX,"
			<< Problem.inputName << ","
			<< Problem.learningRate << ","
			<< elapsed << ","
			<< cplex.getStatus() << endl;
		output.close();
		return;
	}
	elapsed = cplex.getTime() - elapsed;

	ofstream output(Problem.outputName, fstream::app);
	output
		<< "CPLEX,"
		<< Problem.inputName << ","
		<< Problem.learningRate << ","
		<< cplex.getObjValue() << ","
		<< cplex.getBestObjValue() << ","
		<< cplex.getNnodes() << ","
		<< cplex.getNiterations() << ","
		<< cplex.getMIPRelativeGap() << ","
		<< elapsed << endl;
	output.close();
}

SolutionGraph modelo_CP_posicional(Cenario& Problem)
{
	int nop = Problem.nOperations,
		nMachines = Problem.nMachines,
		narc = Problem.nArcs;
	int nIntervalVariables = 0;
	IloEnv env;
	IloModel model(env);

	//declarar variáveis
	IloIntervalVarArray o(env, nop);
	nIntervalVariables += nop;
	IloArray<IloIntervalVarArray2> a(env, nMachines);

	for (IloInt j = 0; j < nMachines; j++) {
		a[j] = IloIntervalVarArray2(env, nop);
		IloInt r = 0;
		for (auto i : Problem.I[j])
		{
			a[j][i] = IloIntervalVarArray(env);
			r++;
		}
	}

	IloNumExprArray ends(env);

	for (IloInt i = 0; i < nop; i++) {
		o[i] = IloIntervalVar(env);
		IloIntervalVarArray members(env);

		for (IloInt k = 0; k < Problem.OperationSet[i].nMachines; k++) {
			IloIntervalVar prec;
			IloInt m = Problem.OperationSet[i].machineSet[k], d = Problem.ProcessingTimes[i][m];

			for (IloInt r = 0; r < (IloInt)Problem.I[m].size(); r++)
			{
				IloIntervalVar member(env, Problem.learningFunction((int)d, (int)r + 1, Problem.learningRate));

				member.setOptional();

				if (prec.getImpl()) {
					model.add(IloEndBeforeStart(env, prec, member));
				}
				prec = member;

				members.add(member);
				a[m][i].add(member);
				nIntervalVariables++;
			}
		}
		model.add(IloAlternative(env, o[i], members));
	}

	for (IloInt k = 0; k < nMachines; k++) {
		for (IloInt r = 0; r < (IloInt)Problem.I[k].size() - 1; r++) {
			IloOr or1(env), or2(env);
			for (auto i : Problem.I[k]) {
				or1.add(IloPresenceOf(env, a[k][i][r]));
				or2.add(IloPresenceOf(env, a[k][i][r + 1]));
			}

			model.add(IloIfThen(env, or2, or1));
		}
	}

	for (IloInt k = 0; k < nMachines; k++) {
		IloIntervalVar prec;
		for (IloInt r = 0; r < (IloInt)Problem.I[k].size() - 1; r++) {
			IloIntervalVarArray members(env);
			for (auto i : Problem.I[k]) {
				for (auto j : Problem.I[k]) {
					model.add(IloEndBeforeStart(env, a[k][i][r], a[k][j][r + 1]));
				}
			}
		}
	}

	for (IloInt a = 0; a < narc; a++) {
		model.add(IloEndBeforeStart(env, o[Problem.PrecedenceArcSet[a].source], o[Problem.PrecedenceArcSet[a].tail]));
		ends.add(IloEndOf(o[Problem.PrecedenceArcSet[a].tail]));
	}

	for (IloInt j = 0; j < nMachines; j++) {
		IloIntervalVarArray members(env);
		for (auto i : Problem.I[j]) {
			members.add(a[j][i]);
		}
		model.add(IloNoOverlap(env, members));
	}

	IloCumulFunctionExpr expr(env);
	for (IloInt i = 0; i < nop; i++) {
		expr += IloPulse(o[i], 1);
	}
	model.add(expr <= nMachines);
	expr.end();

	IloObjective objective = IloMinimize(env, IloMax(ends));
	model.add(objective);

	IloCP cp(model);

	int sum = 0;
	for (int k = 0; k < nMachines; k++)
	{
		sum += (int)Problem.I[k].size();
	}

	cp.setParameter(IloCP::TimeLimit, Problem.timeLimit);
	//cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet);
	//cp.setParameter(IloCP::Workers, 1);
	cp.setParameter(IloCP::RelativeOptimalityTolerance, 0);
	cp.setParameter(IloCP::OptimalityTolerance, 1 - 1e-06);

	bool warm_start = true;

	//TODO
	if (warm_start) {
		SolutionGraph G, Gbest, GSPT, GECT;

		GSPT = SPT(Problem);
		GECT = ECT(Problem);

		G = (GSPT.makespan < GECT.makespan) ? GSPT : GECT;

		IloSolution startsol(env);

		for (IloInt i = 0; i < nop; i++) {
			int i_novo = (int)i + 1;
			startsol.setPresent(a[G.ff[i_novo]][i][G.positions[i_novo] - 1]);
			startsol.setStart(a[G.ff[i_novo]][i][G.positions[i_novo] - 1], G.st[i_novo]);
			startsol.setStart(o[i], (IloInt)G.st[i_novo]);
			startsol.setLength(a[G.ff[i_novo]][i][G.positions[i_novo] - 1], (IloInt)G.weight[i_novo]);
		}

		cp.setStartingPoint(startsol);
	}

	IloNum elapsed;
	elapsed = cp.getTime();

	if (cp.solve()) {
		ofstream instancesInfo("instancesInfo.csv", fstream::app);
		instancesInfo << "CP,"
			<< Problem.inputName << ","
			<< cp.getInfo(IloCP::NumberOfConstraints) << ","
			<< cp.getInfo(IloCP::NumberOfIntervalVariables) << endl;
		instancesInfo.close();
		//return;

		elapsed = cp.getTime() - elapsed;
		cp.out() << "Makespan \t: " << cp.getObjValue() << std::endl;

		ofstream output(Problem.outputName, fstream::app);
		output << "CP,"
			<< Problem.inputName << ","
			<< Problem.learningRate << ","
			<< cp.getObjValue() << ","
			<< cp.getObjBound() << ","
			<< cp.getInfo(IloCP::NumberOfBranches) << ","
			<< cp.getInfo(IloCP::NumberOfConstraints) << ","
			<< cp.getInfo(IloCP::NumberOfIntervalVariables) << ","
			<< cp.getObjGap() << ","
			<< cp.getInfo(IloCP::SolveTime) << ","
			<< elapsed << endl;
		output.close();
	}
	else {
		ofstream output(Problem.outputName, fstream::app);
		output
			<< "CP,"
			<< Problem.inputName << ","
			<< Problem.learningRate << ","
			<< cp.getStatus() << ","
			<< cp.getObjBound() << ","
			<< cp.getInfo(IloCP::NumberOfConstraints) << ","
			<< cp.getInfo(IloCP::NumberOfIntegerVariables) << ","
			<< cp.getInfo(IloCP::SolveTime) << ","
			<< elapsed << endl;
		output.close();
		cp.out() << "No solution found." << std::endl;
	}

	SolutionGraph SolucaoOtima(Problem);

	for (int i = 0; i < SolucaoOtima.nNodes - 2; i++)
	{
		SolucaoOtima.st[i + 1] = (int)cp.getStart(o[i]);
		SolucaoOtima.weight[i + 1] = (int)cp.getSize(o[i]);

		for (auto& k : Problem.F[i])
		{
			for (int r = 0; r < (int)Problem.I[k].size(); r++)
			{
				if (cp.isPresent(a[k][i][r])) {
					SolucaoOtima.ff[i + 1] = k;
					break;
				}
			}
		}
	}
	env.end();

	return SolucaoOtima;
}