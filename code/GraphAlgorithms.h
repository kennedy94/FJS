#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include "Operation.h"

using namespace std;

void CountJobs(vector<bool>& visited, vector<list<int>>& AdjList,
	int vertice, vector<Operation>& op, int job);