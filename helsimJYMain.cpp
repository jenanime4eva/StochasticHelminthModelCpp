/*
Individual-based stochastic simulation of helminths
Created by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/

#include "CSimulator.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <time.h>

using namespace std;

int main(int argc, char** argv)
{
	// Do we have the correct command line?
	if(argc!=4)
	{
		std::cout << "Incorrect command line. Expected 4: Run name(SPACE)Path for log file(SPACE)Path for results file(SPACE)Path for parameter file\n";
		return 0; 
	}

	// Set the seed
	srand(time(NULL)); // Initialise random number generator

	// Create a CSimulator object
	CSimulator simulator;

	// Initialise the input/output aspects of the simulator
	if(!simulator.initialiseIO(argv[1],argv[2],argv[3]))
	{
		simulator.logStream << "Failure in initialiseIO(...).\n" << std::flush;
		return 0;
	}

	simulator.initialiseSimulation(); // General initialisation
	simulator.outputSimulation(); // Output simulation

	return 1; 
}
