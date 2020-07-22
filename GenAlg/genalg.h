#pragma once

#ifdef MAKE_DLL
#define DLL_IMEX __declspec(dllexport)
#else
#define DLL_IMEX __declspec(dllimport)
#endif

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <bitset>
#include <iostream>
#include <thread>

// for dynamic: add DLL_IMEX after keywords 'class' to prepare dynamic dll, and define MAKE_DLL in project settings
extern "C" {
	class DLL_IMEX GeneticAlgorithm {
	private:
		static unsigned __int32 popSize;
		static unsigned __int16 groupSize;
		static unsigned __int16 nGenes;
		static unsigned __int16 generations;
		static unsigned __int32 muRate;
		static unsigned __int32 crRate;

		static bool block;
		static bool multithreaded;

		class DLL_IMEX GaSolution {
		public:
			GaSolution();
			unsigned __int32* genes;
			double fitness;

			void InitGeneString(unsigned __int16 n, bool randomise);
		};

		static GaSolution* population;

	public:
		static GaSolution bestSolution;

		static double(*FitnessEval)(unsigned __int32* inputGenes);

		static void EvaluateFitnessForPop();

		static void CreatePopulation();

		static void GenerateChild(unsigned __int32 posM, unsigned __int32 posD, GaSolution newSol);

		static void TsThread(GaSolution* nextGen, unsigned __int32 start, unsigned __int32 end);

		static void TournamentSelection();

		static void UpdateBest(GaSolution input);

		static void Initialise(__int32 inputPopSize, __int32 inputNumberOfGenerations, __int32 numberOfGenes, double(*fitnessFunction)(unsigned __int32* inGenes), __int32 inputGroupSize, __int32 mutationRateIn100000, __int32 crossoverRateIn100000, bool mt = false);

		static void RunGeneticAlgorithm(bool printOutput);

		static void ClearObject();
	};
}