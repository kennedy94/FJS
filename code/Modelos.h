#pragma once

#include "Cenario.h"
#include "Heuristics.h"
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

//#include "gurobi_c++.h"

void modelo_MILP_posicional_CPLEX(Cenario& Problem);

SolutionGraph modelo_CP_posicional(Cenario& Problem);