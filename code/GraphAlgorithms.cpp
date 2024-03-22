#include "GraphAlgorithms.h"

void CountJobs(vector<bool>& visited, vector<list<int>>& AdjList, int v, vector<Operation>& operation, int idJob) {
	operation[v].idJob = idJob;
	visited[v] = true;
	for (auto& adj : AdjList[v])
		if (!visited[adj])
			CountJobs(visited, AdjList, adj, operation, idJob);
}