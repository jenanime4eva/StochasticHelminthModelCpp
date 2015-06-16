/*
Individual-based stochastic simulation of helminths
Created by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/

#include "CSimulator.h"
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
		std::cout << "Incorrect command line. 4 expected arguments: .exe file(space)Run name e.g. Test(space)Path to output log and results files(space)Path of parameter file\n";
		return 0; 
	}

	int t1 = time(NULL); // Start clock counter

	// Set the seed
	// Initialise random number generator
	//srand(time(NULL));
	srand(1); // Fixed seed for testing

	// Create a CSimulator object
	CSimulator simulator;

	// Initialise the input/output aspects of the simulator
	if(!simulator.initialiseIO(argv[1],argv[2],argv[3]))
	{
		simulator.logStream << "Failure in initialiseIO(...).\n" << std::flush;
		return 0;
	}

	simulator.initialiseSimulation(); // General initialisation
	simulator.runSimulation(); // Run simulation
	simulator.outputSimulation(); // Output simulation

	int t2 = time(NULL); // Stop clock counter

	printf ("\nTime taken = %d secs\n", t2 - t1); // Print time taken in seconds to run the program to the console

	return 1; 
}
