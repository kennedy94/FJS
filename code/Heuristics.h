#pragma once
#include "SolutionGraph.h"
#include "GraphAlgorithms.h"
#include <set>
#include <iomanip>
#include <sstream>

const bool IMPRIMIR_ARQUIVO = true;

SolutionGraph SPT(Cenario& FJSSFLE);

SolutionGraph ECT(Cenario& FJSSFLE);

SolutionGraph ReducedLocalSearch(Cenario& FJSSFLE, SolutionGraph& G, string localSearchStrategy, std::chrono::high_resolution_clock::time_point tic, bool onlyCriticalOperations = false);

SolutionGraph LocalSearch(Cenario& FJSSFLE, SolutionGraph& G, string localSearchStrategy, double tolerancy, std::chrono::high_resolution_clock::time_point tic, bool onlyCriticalOperations = false);

SolutionGraph ILS(Cenario& FJSSFLE, string localSearchStrategy, int pertMin, int pertMax, double tolerancy, bool tuning = false, bool onlyCriticalOperations = false);

SolutionGraph SA(Cenario& FJSSFLE, string localSearch, string localSearchStrategy, int itMin, int itMax, double tempM, double tempP, double tempF, double deltaMin, double deltaMax, double tolerancy, string functionType, bool tuning = false);

SolutionGraph TS(Cenario& FJSSFLE, int listSizeMax, bool tuning = false, bool onlyCriticalOperations = false);

SolutionGraph GRASP(Cenario& FJSSFLE, string localSearchStrategy, double alphaGrasp, double tolerancy, bool tuning = false, bool onlyCriticalOperations = false);

SolutionGraph RandomizedECT(Cenario& FJSSFLE, double alphaGrasp);

SolutionGraph RandomizedSPT(Cenario& FJSSFLE, double alphaGrasp);
