#pragma once

class SolutionNode
{
public:
	int id;
	int label;
	double weight;

	SolutionNode(int id, int label, double weight);

	bool operator ==(SolutionNode B);
};
