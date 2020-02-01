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
	void GeneticAlgorithm::GASolution::InitGeneString(unsigned __int16 n, bool randomise) {
		genes = new unsigned __int32[n];

		if (randomise) {
			srand(time(NULL));

			for (unsigned __int16 i = 0; i < n; i++) {
				genes[i] = unsigned __int32(rand() % 4294967296);
			}
		}
	}

	GeneticAlgorithm::GeneticAlgorithm(__int32 inputPopSize, __int32 inputNumberOfGenerations, __int32 numberOfGenes, double(*fitnessFunction)(unsigned __int32* inGenes), __int32 inputGroupSize, __int32 mutationRateIn100000, __int32 crossoverRateIn100000) {
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

		for (unsigned __int32 i = 0; i < popSize; i++) {
			population[i].InitGeneString(nGenes, true);
		}

		bestSolution.InitGeneString(nGenes, false);
		bestSolution.fitness = -50000000.0;
	}

	void GeneticAlgorithm::EvaluateFitnessForPop() {
		for (unsigned __int32 i = 0; i < popSize; i++) {
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

	void GeneticAlgorithm::GenerateChild(unsigned __int32 posM, unsigned __int32 posD, GeneticAlgorithm::GASolution newSol) {
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

		for (unsigned __int32 i = 0; i < nGenes; i++) {
			newSol.genes[i] = std::bitset<32>(childbin.substr(i * 32, 32)).to_ulong();
		}
	}

	void GeneticAlgorithm::TournamentSelection() {
		GASolution* nextGen = new GASolution[popSize];
		for (unsigned __int32 i = 0; i < popSize; i++) {
			nextGen[i].InitGeneString(nGenes, false);
		}

		bool unique;
		unsigned __int32 bestContender;
		unsigned __int32* contenders = new unsigned __int32[groupSize];


		for (unsigned __int32 popMember = 0; popMember < popSize; popMember++) {

			// for m
			for (int i = 0; i < groupSize; i++) {
				contenders[i] = 0;
			}

			for (int i = 0; i < groupSize; i++) {
				unique = false;
				while (!unique) {
					contenders[i] = unsigned __int32(rand() % popSize);
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
					contenders[i] = unsigned __int32(rand() % popSize);
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

			GenerateChild(mother, father, nextGen[popMember]);
		}

		for (unsigned __int32 i = 0; i < popSize; i++) {
			for (unsigned __int32 j = 0; j < nGenes; j++) {
				population[i].genes[j] = nextGen[i].genes[j];				
			}
			delete[] nextGen[i].genes;
		}
		delete[] nextGen;
	}

	void GeneticAlgorithm::RunGeneticAlgorithm(bool printOutput) {
		if (block) {
			std::cout << "\nGenAlg instance not initialised properly, use the parametised constructor only!\n";
		}
		CreatePopulation();

		std::cout.precision(15);
		for (int i = 0; i < generations; i++) {
			EvaluateFitnessForPop();
			if (printOutput) std::cout << std::scientific << "Generation " << i + 1 << " completed. Best fitness = " << bestSolution.fitness << std::endl;
			TournamentSelection();
		}
	}

	void GeneticAlgorithm::ClearObject() {
		for (unsigned __int32 i = 0; i < popSize; i++) {
			delete[] population[i].genes;
		}
		delete[] population;
	}

}