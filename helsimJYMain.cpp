/*
Individual-based stochastic simulation of helminths
Written by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include "randlib.h"
#include "SimulatorJY.h"  // Include the CSimulator class

int main(int argc, char** argv)
{
	// Do we have the correct command line?
	if(argc!=4)
	{
		std::cout << "Incorrect command line. Expected:\nFile.exe logFile paramFile resultsFile\n";
		return 0; 
	}

	// Set the seed
	setall(1,1); 

	// Create a CSimulator object
	helsimJY mySim;

	mySim.initialiseIO(argv[1], argv[2], argv[3]); 
	
	mySim.initialiseSimulation(); 

	mySim.runRealizations(); 

	mySim.outputRealisation(0);

	return 1; 
}
