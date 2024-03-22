#include <iostream>
#include <chrono>
#include <random>
#include <string>
#include <fstream>
#include <string.h>
#include "Cenario.h"
#include "Heuristics.h"
#include "Modelos.h"
#include "FJS.h"
#include "Tayebi.h"
#include <iomanip>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 5) {
		cerr << "Erro de input" << endl;
		exit(0);
	}

	//Read Parameters
	const char
		* inputName = argv[1 * 2],
		* outputName = argv[2 * 2];
	const double
		timeLimit = stod(argv[3 * 2]),
		learningRate = stod(argv[4 * 2]);
	int SEED = stoi(argv[5 * 2]);
	if (SEED == 0)
		SEED = (int)time(NULL);

#ifdef _DEBUG
	cout << "Instance: " << inputName << endl << "Output: " << outputName << endl;
	cout << "Time limit: " << timeLimit << endl << "Learning rate: " << learningRate << endl;
#endif
	bool tuning = false;
	if (strcmp(outputName, "tuning.out") == 0)
		tuning = true;

	/*Cenario* Teste = new Cenario(inputName, outputName, learningRate, timeLimit, SEED);
	return -1;*/
	SolutionGraph G;

	if (strcmp(argv[6 * 2], "MODELOMILP") == 0)
	{
#ifdef _DEBUG
		cout << "Model: " << argv[6 * 2] << endl;
#endif
		Cenario FJSSFLE_MODELO(inputName, outputName, learningRate, timeLimit, SEED);
		modelo_MILP_posicional_CPLEX(FJSSFLE_MODELO);
		return 0;
	}

	if (strcmp(argv[6 * 2], "MODELOCP") == 0)
	{
#ifdef _DEBUG
		cout << "Model: " << argv[6 * 2] << endl;
#endif
		Cenario FJSSFLE_MODELO(inputName, outputName, learningRate, timeLimit, SEED);
		//return -1;
		G = modelo_CP_posicional(FJSSFLE_MODELO);
#ifdef _DEBUG
		G.Feasible(FJSSFLE_MODELO);
		stringstream stream;
		stream << fixed << setprecision(2) << learningRate;
		string method = stream.str();

		vector<int> jobs(FJSSFLE_MODELO.nOperations, 0);
		int cont = 1;
		for (auto& i : FJSSFLE_MODELO.Jobs) {
			for (auto& j : i) {
				jobs[j] = cont;
			}
			cont++;
		}
		G.printGanttChart(inputName, method, "ModeloCP", jobs);
#endif
		return 0;
	}

	int makespan = INT_MAX;
	if (strcmp(argv[6 * 2], "TAYEBI") == 0) {
#ifdef _DEBUG
		cout << "Metaheuristic: " << argv[6 * 2] << endl;
#endif
		FJS FJSSFLE_MODELO(inputName, outputName, timeLimit, learningRate);
		auto tic = std::chrono::high_resolution_clock::now();
		auto solution = TayebiGA(FJSSFLE_MODELO, 60, 3, 0.3, 0.2, 0.5, 6);
		makespan = solution.fitness;
		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;
		ofstream output_ofstream;
		output_ofstream.open(outputName, iostream::app);
		output_ofstream
			<< argv[6 * 2] << ","
			<< inputName << ","
			<< learningRate << ","
			<< makespan << ","
			<< toc << ","
			<< solution.time_to_find << endl;
		output_ofstream.close();

		return 0;
	}

	const char
		* localSearch = argv[6 * 2],			// Full - Reduced - CriticalReduced
		* localSearchStrategy = argv[7 * 2],	// Best - First - Alternate
		* constructiveHeuristic = argv[8 * 2],	// SPT - ECT - Best
		* metaHeuristic = argv[9 * 2];			// ILS - SAlog - SAlin - SAsig - TS - GRASP

#ifdef _DEBUG
	cout << "Local search: " << localSearch << endl << "Local search strategy: " << localSearchStrategy << endl;
	cout << "Constructive heuristic: " << constructiveHeuristic << endl << "Metaheuristic: " << metaHeuristic << endl;
#endif // _DEBUG

	int itMax = 0, pertMin = 0, pertMax = 0, listSizeMax = 0;
	double tempM = 0.0, tempP = 0.0, tempF = 0.0, alphaGrasp = 0.0, deltaMin = 0.0, deltaMax = 0.0, tolerancy = 0.0;
	bool onlyCriticalOperations = false;
	string functionType = "fix";

	if (strcmp(localSearch, "None") != 0) {
		tolerancy = stod(argv[10 * 2]);
		onlyCriticalOperations = stoi(argv[11 * 2]);

#ifdef _DEBUG
		cout << "Tolerancy: " << tolerancy << endl << "Critical operations ? " << onlyCriticalOperations << endl;
#endif // _DEBUG
	}

	if (strcmp(metaHeuristic, "None") != 0)
		itMax = stoi(argv[12 * 2]);

	if (strcmp(metaHeuristic, "ILS") == 0) {
		pertMin = stoi(argv[13 * 2]);
		pertMax = stoi(argv[14 * 2]);
	}
	if (strcmp(metaHeuristic, "SA") == 0) {
		pertMin = stoi(argv[13 * 2]) * 1000;
		pertMax = stoi(argv[14 * 2]) * 1000;
		tempM = stod(argv[15 * 2]);
		tempP = stod(argv[16 * 2]);
		tempF = stod(argv[17 * 2]);
		deltaMin = stod(argv[18 * 2]);
		deltaMax = stod(argv[19 * 2]);
		if (0 > stoi(argv[14 * 2]))
			pertMax = pertMin;
		if (0 > stoi(argv[19 * 2]))
			deltaMax = deltaMin;
		functionType = argv[20 * 2];
	}
	if (strcmp(metaHeuristic, "TS") == 0) {
		listSizeMax = stoi(argv[13 * 2]);
	}
	if (strcmp(metaHeuristic, "GRASP") == 0) {
		alphaGrasp = stod(argv[13 * 2]);
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Cenario* FJSSFLE;
	try
	{
		FJSSFLE = new Cenario(inputName, outputName, learningRate, timeLimit, SEED);
		//return -1;
	}
	catch (const std::exception&)
	{
		cerr << "Erro leitura de instância" << endl;
		ofstream output_ofstream;
		output_ofstream.open("ErroLog.txt", iostream::app);
		output_ofstream << inputName << "," << SEED << "," << metaHeuristic << "," << itMax << "," << pertMin << "," << pertMax << "," << learningRate << endl;
		output_ofstream.close();
		exit(-1);
	}

	std::chrono::high_resolution_clock::time_point tic = std::chrono::high_resolution_clock::now();

	if (strcmp(metaHeuristic, "None") == 0) {
		if (strcmp(constructiveHeuristic, "SPT") == 0)
			G = SPT(*FJSSFLE);
		if (strcmp(constructiveHeuristic, "ECT") == 0)
			G = ECT(*FJSSFLE);
		if (strcmp(constructiveHeuristic, "Best") == 0) {
			SolutionGraph GSPT = SPT(*FJSSFLE);
			SolutionGraph GECT = ECT(*FJSSFLE);
			G = (GSPT.makespan < GECT.makespan) ? GSPT : GECT;
		}

		if (strcmp(localSearch, "None") != 0) {
			G = LocalSearch(*FJSSFLE, G, localSearchStrategy, tolerancy, tic, onlyCriticalOperations);
		}
	}

	try
	{
		if (strcmp(metaHeuristic, "ILS") == 0) {
			ofstream output_ofstream;
			output_ofstream.open("RunningLog.txt", iostream::app);
			output_ofstream << inputName << "," << SEED << "," << metaHeuristic << "," << itMax << "," << pertMin << "," << pertMax << "," << learningRate << endl;
			output_ofstream.close();

			G = ILS(*FJSSFLE, localSearchStrategy, pertMin, pertMax, tolerancy, tuning, onlyCriticalOperations);
		}
		if (strcmp(metaHeuristic, "SA") == 0) {
			G = SA(*FJSSFLE, localSearch, localSearchStrategy, pertMin, pertMax, tempM, tempP, tempF, deltaMin, deltaMax, tolerancy, functionType, tuning);
		}
		if (strcmp(metaHeuristic, "TS") == 0) {
			G = TS(*FJSSFLE, listSizeMax, tuning, onlyCriticalOperations);
		}
		if (strcmp(metaHeuristic, "GRASP") == 0) {
			G = GRASP(*FJSSFLE, localSearchStrategy, alphaGrasp, tolerancy, tuning, onlyCriticalOperations);
		}

		auto toc = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() / 1000.0;

		//system("mkdir output");
		ofstream output_ofstream;
		output_ofstream.open(outputName, iostream::app);
		if (!G.Feasible(*FJSSFLE)) {
			output_ofstream << metaHeuristic << "," << localSearch << "," << localSearchStrategy << "," << constructiveHeuristic << ","
				<< inputName << "," << learningRate << "," << "Erro Inviável" << "," << toc << "," << G.timeToGetSolution << endl;
			output_ofstream.close();
			return 0;
		}

		if (strcmp(metaHeuristic, "None") == 0) {
			output_ofstream << metaHeuristic << "," << localSearch << "," << localSearchStrategy << "," << constructiveHeuristic << ","
				<< inputName << "," << learningRate << "," << G.makespan << "," << toc << "," << G.timeToGetSolution << "," << G.LocalSearchIterations << "," << G.nNeighbours << endl;
		}

		if (strcmp(metaHeuristic, "ILS") == 0) {
			output_ofstream << metaHeuristic << "," << localSearch << "," << localSearchStrategy << "," << constructiveHeuristic << "," << itMax << "," << pertMin << "," << pertMax << ","
				<< inputName << ", " << learningRate << ", " << G.makespan << ", " << toc << "," << G.timeToGetSolution << endl;
		}
		if (strcmp(metaHeuristic, "SA") == 0) {
			output_ofstream << metaHeuristic << "," << localSearch << "," << localSearchStrategy << "," << constructiveHeuristic << "," << pertMin << "," << pertMax << "," << tempM << ","
				<< tempP << "," << tempF << "," << deltaMin << "," << deltaMax << "," << inputName << ", " << learningRate << ", " << G.makespan << "," << toc << "," << G.timeToGetSolution << endl;
		}
		if (strcmp(metaHeuristic, "TS") == 0) {
			output_ofstream << metaHeuristic << "," << localSearch << "," << localSearchStrategy << "," << constructiveHeuristic << "," << itMax << "," << listSizeMax << ","
				<< inputName << ", " << learningRate << ", " << G.makespan << ", " << toc << "," << G.timeToGetSolution << endl;
		}
		if (strcmp(metaHeuristic, "GRASP") == 0) {
			output_ofstream << metaHeuristic << "," << localSearch << "," << localSearchStrategy << "," << constructiveHeuristic << "," << itMax << "," << alphaGrasp << ","
				<< inputName << ", " << learningRate << ", " << G.makespan << ", " << toc << "," << G.timeToGetSolution << endl;
		}

		output_ofstream.close();

		//cout << "cost: " << G.makespan << endl;
		cout << G.makespan;
		//cout << "time: " << toc << endl;
	}
	catch (const std::exception& e)
	{
		cout << e.what() << endl;
		ofstream output_ofstream;
		output_ofstream.open("ErrorFile.err", iostream::app);

		output_ofstream << "Erro" << e.what() << "," << SEED << "," << metaHeuristic << "," << itMax << "," << pertMin << "," << pertMax << "," << inputName << ", " << learningRate << "," << G.makespan << endl;

		output_ofstream.close();
	}

#ifdef _DEBUG
	if (G.Feasible(*FJSSFLE))
		cout << "\n\nViavel!" << endl;
	else {
		cout << "\n\nInviavel" << endl;
		exit(-1);
	}
	//return 0;
	//G.UpdateStartingTimes(FJSSFLE);
	stringstream stream;
	stream << fixed << setprecision(2) << learningRate;
	string method = stream.str();

	vector<int> jobs(FJSSFLE->nOperations, 0);
	int cont = 1;
	for (auto& i : FJSSFLE->Jobs) {
		for (auto& j : i) {
			jobs[j] = cont;
		}
		cont++;
	}
	G.printGanttChart(inputName, method, constructiveHeuristic, jobs);
#endif

	return 0;
}