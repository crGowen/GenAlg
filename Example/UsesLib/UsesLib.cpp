// UsesLib.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <genalg.h>

double EvalSolut(unsigned __int32* g) {

	double x = (double(g[0]) / 4294967295.0) * 20.0;
	double y = (double(g[1]) / 4294967295.0) * 20.0;
	double res = 15.0 - x * y;

	res = res * res;

	return res*(-1.0);
}

int main()
{
	GenAlg::GeneticAlgorithm ga(100000, 100, 2, &EvalSolut, 2, 10000, 40000);
	ga.RunGeneticAlgorithm();

	for (int i = 0; i < 2; i++) {
		std::cout << "Gene " << i << " = " << (double(ga.bestSolution.genes[i]) / 4294967295.0) * 20.0 << "\n";
	}
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
