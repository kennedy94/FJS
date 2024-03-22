#include "SolutionNode.h"

SolutionNode::SolutionNode(int id, int label, double weight) {
	this->id = id;
	this->label = label;
	this->weight = weight;
}

bool SolutionNode::operator==(SolutionNode B)
{
	return this->id == B.id;
}