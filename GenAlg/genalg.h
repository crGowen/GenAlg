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

			~GASolution();

			void InitGeneString(unsigned __int16 n, bool randomise);
		};
	public:
		GASolution* population;
		GASolution bestSolution;

		GeneticAlgorithm(unsigned __int32 inputPopSize, unsigned __int16 inputNumberOfGenerations, unsigned __int16 numberOfGenes, double(*fitnessFunction)(unsigned __int32* inGenes), unsigned __int16 inputGroupSize, unsigned __int32 mutationRateIn100000, unsigned __int32 crossoverRateIn100000);

		double(*FitnessEval)(unsigned __int32* inputGenes);

		void EvaluateFitnessForPop();

		void CreatePopulation();

		GASolution GenerateChild(unsigned __int32 posM, unsigned __int32 posD);
		
		void TournamentSelection();

		void UpdateBest(GASolution input);

		void RunGeneticAlgorithm();
	};
}