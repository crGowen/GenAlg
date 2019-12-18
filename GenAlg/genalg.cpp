#include "stdafx.h"
#include "genalg.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <bitset>
#include <iostream>

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

	GeneticAlgorithm::GeneticAlgorithm(unsigned __int32 inputPopSize, unsigned __int16 inputNumberOfGenerations, unsigned __int16 numberOfGenes, double(*fitnessFunction)(unsigned __int32* inGenes), unsigned __int16 inputGroupSize, unsigned __int32 mutationRateIn100000, unsigned __int32 crossoverRateIn100000) {
		popSize = inputPopSize;
		generations = inputNumberOfGenerations;
		nGenes = numberOfGenes;
		FitnessEval = fitnessFunction;
		groupSize = inputGroupSize;			
		muRate = mutationRateIn100000;
		crRate = crossoverRateIn100000;

		block = false;
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
		std::string mbin(32 * nGenes, '2');
		for (unsigned __int32 i = 0; i < nGenes; i++) {
			mbin.replace(i * 32, 32, std::bitset<32>(population[posM].genes[i]).to_string());
		}

		std::string dbin(32 * nGenes, '2');
		for (unsigned __int32 i = 0; i < nGenes; i++) {
			dbin.replace(i * 32, 32, std::bitset<32>(population[posD].genes[i]).to_string());
		}

		unsigned __int32 tempRand;
		std::string childbin = mbin;

		if (rand() % 100000 < crRate) {
			//crossover
			tempRand = (rand() % (childbin.length() - 1)) + 1;
			childbin.replace(tempRand, childbin.length() - tempRand, dbin.substr(tempRand, childbin.length() - tempRand));
		}

		if (rand() % 100000 < muRate) {
			tempRand = rand() % childbin.length();
			if (childbin[tempRand] == '0') childbin.replace(tempRand, 1, "1");
			else childbin.replace(tempRand, 1, "0");
		}

		GeneticAlgorithm::GASolution output;
		for (unsigned __int32 i = 0; i < nGenes; i++) {
			output.genes[i] = std::bitset<32>(childbin.substr(i * 32, 32)).to_ulong();
		}

		return output;
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

		for (unsigned __int32 i = 0; i < popSize; i++) {
			for (unsigned __int32 j = 0; j < nGenes; i++) {
				population[i].genes[j] = nextGen[i].genes[j];				
			}
			delete[] nextGen[i].genes;
		}
		delete[] nextGen;
	}

	void GeneticAlgorithm::RunGeneticAlgorithm() {
		CreatePopulation();

		for (int i = 0; i < generations; i++) {
			EvaluateFitnessForPop();
			std::cout << "Generation " << i + 1 << " completed. Best fitness = " << bestSolution.fitness << std::endl;
			TournamentSelection();
		}
	}

}