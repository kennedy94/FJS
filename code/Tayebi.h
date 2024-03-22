#pragma once
#include "FJS.h"
#include<iostream>
#include <random>
#include <set>
#include <vector>
#include <algorithm>

#include <chrono>

struct Chromosome {
	vector< int>
		operations,
		machines;
	int fitness = 100000000;
	double time_to_find = 0.0;

	bool operator<(Chromosome a) {
		return this->fitness < a.fitness;
	}

	bool operator==(Chromosome a) {
		for (int i = 0; i < (int)a.machines.size(); i++) {
			if (a.machines[i] != this->machines[i] || a.operations[i] != this->operations[i])
				return false;
		}
		return true;
	}
};

Chromosome TayebiGA(FJS& P, int npop, int M,
	double Pc, double Pns, double Paf, int Lmax);

void evaluateChromosome(FJS& P, Chromosome& S);

void repairChromosome(FJS& P, Chromosome& S);

vector<Chromosome> generatePopulation(FJS& P, int popsize);

Chromosome neighbourhood(FJS& P, Chromosome& S, int movimento);

Chromosome insertion(FJS& P, Chromosome& S, string scenario, int i, int j);

Chromosome swap(FJS& P, Chromosome& S, string scenario, int i, int j);

Chromosome inversion(FJS& P, Chromosome& S, string scenario, int i, int j);

Chromosome threepoint(FJS& P, Chromosome& S, string scenario);

Chromosome fourpoint(FJS& P, Chromosome& S, string scenario);

vector<Chromosome> xSelections(vector<Chromosome>& popu, int M);

void repairChromosome(FJS& P, Chromosome& S);

void crossover(FJS& P, vector<Chromosome>& pop2,
	Chromosome& P1, Chromosome& P2,
	string mode, string scenario);

bool isFeasible(FJS& P, Chromosome S);