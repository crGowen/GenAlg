namespace GenAlg {
	class GeneticAlgorithm {
	private:
		unsigned __int32 popSize;
		unsigned __int16 groupSize;
		unsigned __int16 nGenes;
		unsigned __int16 generations;
		unsigned __int32 muRate;
		unsigned __int32 crRate;

		bool block = true;

		class GASolution {
		public:
			unsigned __int32* genes;
			double fitness;

			void InitGeneString(unsigned __int16 n, bool randomise);
		};

		GASolution* population;

	public:
		GASolution bestSolution;

		double(*FitnessEval)(unsigned __int32* inputGenes);

		void EvaluateFitnessForPop();

		void CreatePopulation();

		void GenerateChild(unsigned __int32 posM, unsigned __int32 posD, GASolution newSol);

		void TournamentSelection();

		void UpdateBest(GASolution input);

		GeneticAlgorithm(__int32 inputPopSize, __int32 inputNumberOfGenerations, __int32 numberOfGenes, double(*fitnessFunction)(unsigned __int32* inGenes), __int32 inputGroupSize, __int32 mutationRateIn100000, __int32 crossoverRateIn100000);

		void RunGeneticAlgorithm(bool printOutput);

		void ClearObject();
	};
}