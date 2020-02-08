// UsesLib.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <genalg.h>

// THIS PROGRAM WILL ONLY COMPILE FOR x64/Release unless you also recompile the GenAlg library for other platforms
// ALSO NOTE: the algorithm is not necessarily going to produce the most 'elegant' solution, it will simply go for a solution that works well, for instance: in this example program,
// the obvious solution to this fitness function (produce 15 from x * y) is x=3, y=5, however a result of something like
// x = 2.6660446, y = 5.6263123 is just as likely as they both produce 15 within the floating point margin of error

// fitness function must always be defined by the user as a double with a (unsigned __int32*) parameter
double EvalSolut(unsigned __int32* g) {
	double x = (double(g[0]) / 4294967295.0) * 10.0; // convert to double and divide by 4294967295.0 to get a range of ~0-2 (also with the negatives: -2 to 2), then multiply by 10 to get a range of 0-20 (-20 to 20) with high 'resolution'
	double y = (double(g[1]) / 4294967295.0) * 10.0;

	// above two lines get the genes and translate them into numbers ranging from -20 to 20

	double res = 15.0 - x * y; // our fitness function wants to use the two genes as input to get as close to 15 as possible

	if (res>0) res *= (-1.0);
	// apply a negative sign (because we want to minimise the difference from 15, and the fitness values are interpreted as 'greater is better',
	// so by making them negative, a 'smaller' -0.000002 is 'greater' than '-2' because they are negative

	return res;
}

int main()
{
	// create algorithm with:
	/*
	a population of 60,000
	200 generations
	2 genes (which are x and y in our fitness evaluation that we defined above)
	we point to our defined fitness evaluation "EvalSolut",
	tournament selection group size of 2
	20,000/100,000 (20%) mutation rate
	60,000/100,000 (60%) crossover rate
	*/
	GenAlg::GeneticAlgorithm ga(60000, 200, 2, &EvalSolut, 2, 20000, 60000);

	//run the algorithm with output as true (output means it will print like 'Generation 80 completed. Best fitness = -4.21e-19')
	ga.RunGeneticAlgorithm(true);

	// now that alg is finished, grab best solution and print the value of its genes (use the same translation as in the fitness eval for consistency)
	for (int i = 0; i < 2; i++) {
		std::cout << "Gene " << i << " = " << (double(ga.bestSolution.genes[i]) / 4294967295.0) * 10.0 << "\n";
	}

	// clear dynamic memory used
	ga.ClearObject();
}
