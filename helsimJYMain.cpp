/*
Individual-based stochastic simulation of helminths
Created by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/

//#include "CParamReader.h"
#include ".\CSimulator.h"
//#include "CHost.h"
#include "CPreDetEventQueue.h"

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
		std::cout << "Incorrect number of command line arguments. Expected 4: run name, path for results/log and param file with path\n";
		return 0; 
	}


	// Set the seed.
	//srand(10);
	srand(time(NULL));

	// Create a CSimulator object
	CSimulator simulator;
	// Initialise the input/output aspects of the simulator

	if(!simulator.initialiseIO(argv[1],argv[2],argv[3]))
	{
		simulator.logStream << "Failure in initialiseIO(...).\n" << std::flush;
		return 0;
	}


	simulator.initialiseSimulation(); // General initialisation
	//simulator.runSimulation(); // Run simulation
	simulator.outputSimulation(); // Output simulation

	return 1; 
}
