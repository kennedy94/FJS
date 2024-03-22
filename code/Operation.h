#pragma once
#include <vector>
class Operation
{
public:
	int
		idOperation,
		idJob,
		nMachines;
	std::vector< int> machineSet;
	std::vector<int> processingTimes;

	Operation() {
		idOperation = 0;
		idJob = 0;
		nMachines = 0;
	}
};
