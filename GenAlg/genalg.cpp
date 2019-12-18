#include "stdafx.h"
#include "genalg.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

namespace GenAlg
{
	GeneticAlgorithm::GASolution::~GASolution() {
		delete[] genes;
	}

	void GeneticAlgorithm::GASolution::InitGeneString(unsigned __int16 n, bool randomise) {
		genes = new unsigned __int32[n];

		if (randomise) {
			srand(time(NULL));

			for (unsigned __int16 i = 0; i < n; i++) {
				genes[i] = unsigned __int32(rand() % 4294967296);
			}
		}
	}

	void GeneticAlgorithm::CreatePopulation() {
		population = new GASolution[popSize];

		for (unsigned __int32 i; i < popSize; i++) {
			population[i].InitGeneString(nGenes, true);
		}

		bestSolution.InitGeneString(nGenes, false);
	}

	void GeneticAlgorithm::EvaluateFitnessForPop() {
		for (unsigned __int32 i; i < popSize; i++) {
			population[i].fitness = FitnessEval(population[i].genes);
			if (population[i].fitness > bestSolution.fitness) {
				UpdateBest(population[i]);
			}
		}
	}

	void GeneticAlgorithm::UpdateBest(GeneticAlgorithm::GASolution input) {
		bestSolution.fitness = input.fitness;
		for (int i = 0; i < nGenes; i++) {
			bestSolution.genes[i] = input.genes[i];
		}
	}

	GeneticAlgorithm::GASolution GeneticAlgorithm::GenerateChild(unsigned __int32 posM, unsigned __int32 posD) {
		//DO
	}

	void GeneticAlgorithm::TournamentSelection() {
		GASolution* nextGen = new GASolution[popSize];
		bool unique;
		unsigned __int32 bestContender;
		unsigned __int32* contenders = new unsigned __int32[groupSize];


		for (unsigned __int32 popMember; popMember < popSize; popMember++) {

			// for m
			for (int i = 0; i < groupSize; i++) {
				contenders[i] = 0;
			}

			for (int i = 0; i < groupSize; i++) {
				unique = false;
				while (!unique) {
					contenders[i] = unsigned __int32(rand() % 4294967296);
					unique = true;
					for (int j = 0; j < groupSize; j++) {
						if (contenders[i] == contenders[j] && i != j) {
							unique = false;
						}
					}
				}
			}
			//initialised contenders, now find the best 1
			bestContender = contenders[0];
			for (int i = 1; i < groupSize; i++) {
				if (population[contenders[i]].fitness > population[bestContender].fitness) {
					bestContender = contenders[i];
				}
			}

			unsigned __int32 mother = bestContender;

			// for d
			for (int i = 0; i < groupSize; i++) {
				contenders[i] = 0;
			}

			for (int i = 0; i < groupSize; i++) {
				unique = false;
				while (!unique) {
					contenders[i] = unsigned __int32(rand() % 4294967296);
					unique = !(contenders[i] == mother);
					for (int j = 0; j < groupSize; j++) {
						if (contenders[i] == contenders[j] && i != j) {
							unique = false;
						}
					}
				}
			}

			//initialised contenders, now find the best 1
			bestContender = contenders[0];
			for (int i = 1; i < groupSize; i++) {
				if (population[contenders[i]].fitness > population[bestContender].fitness) {
					bestContender = contenders[i];
				}
			}

			unsigned __int32 father = bestContender;

			nextGen[popMember] = GenerateChild(mother, father);
		}

		//now set current generation to nextgen
	}

	void GeneticAlgorithm::RunGeneticAlgorithm(unsigned __int32 inputPopSize, unsigned __int16 inputGens, unsigned __int16 numberOfGenes, double(*fitnessFunction)(unsigned __int32* inGenes), unsigned __int16 inGS) {
		popSize = inputPopSize;
		groupSize = inGS;
		nGenes = numberOfGenes;
		generations = inputGens;
		FitnessEval = fitnessFunction;

		CreatePopulation();

		for (int i = 0; i < generations; i++) {
			EvaluateFitnessForPop();
			// Say "Gen complete" with best solution
			TournamentSelection();
		}
	}

}