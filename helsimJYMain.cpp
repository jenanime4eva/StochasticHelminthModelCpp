/*
Individual-based stochastic simulation of helminths
Created by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/

//#include "CParamReader.h"
#include ".\CSimulator.h"
//#include "CHost.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <random>

using namespace std;

int main(int argc, char** argv)
{
	// Do we have the correct command line?
	if(argc!=4)
	{
		std::cout << "Incorrect command line. Expected arguments for your .exe file: logFile paramFile resultsFile\n";
		return 0; 
	}

	// Set the seed (ADD THIS LATER)

	// Create a CSimulator object
	CSimulator SIMULATE;
	SIMULATE.initialiseIO(argv[1],argv[2],argv[3]); // Initialise the input/output aspects of the simulator
	SIMULATE.initialiseSimulation(); // General initialisation
	SIMULATE.runSimulation(); // Run simulation
	SIMULATE.outputSimulation(0); // Output simulation

	return 1; 
}
