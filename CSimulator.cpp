/*
 * CSimulator.cpp
 *
 *  Created on: 26 Nov 2014
 *      Author: Jie Yang
 */

#include ".\CSimulator.h"
#include "CHost.h"
#include "CParamReader.h"
#include <iostream>
#include <string.h>
#include <stdlib.h>


// Class Constructor
CSimulator::CSimulator()
{
	hostPopulation = NULL;
	results = NULL;
	nRepetitions = 0;
}

// Class Destructor
CSimulator::~CSimulator()
{
	// Delete results allocations
	if (results!=NULL)
	{
		// Delete memory
		for (int i=0;i<nRepetitions;i++)
		{
			delete results[i];
		}

		delete[] results;
	}

	// Delete hostPopulation array.
	if (hostPopulation!=NULL)
	{
		for (int i=0;i<nHosts;i++)
		{
			delete hostPopulation[i];
		}
		delete[] hostPopulation;
	}
}

// Initialise the input/output aspects of the simulator
bool CSimulator::initialiseIO(char* logFileName, char* paramFileName, char* resultsFileName)
{
	// Setting up all the files needed
	logStream.open(logFileName);
	if(!logStream.is_open())
	{
		std::cout << "Couldn't open the log file: " << logFileName << "\nexiting\n";
		return false;
	}
	logStream << "Log file opened.\n";

	if (!myReader.setNewFileName(paramFileName))
 	{
		logStream << "Couldn't open the param file: " << paramFileName << "\nexiting\n";
		return false;
	}

	resultsStream.open(resultsFileName);
	if(!resultsStream.is_open())
	{
		logStream << "Couldn't open the results file: " << resultsFileName << "\nexiting\n";
		return false;
	}


	nRepetitions = atoi(myReader.getParamString("param1"));
	//nYears = atoi(myReader.getParamString("param2"));
	//nOutputsPerYear = atoi(myReader.getParamString("param3"));

	// Everything is ok
	return true;
}

// General initialisation
bool CSimulator::initialiseSimulation()
{
	// Get parameters
	nHosts = atoi(myReader.getParamString("param4"));

	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;
	}

	return true;
}

// Run simulation
void CSimulator::runSimulation()
{

}

// Output simulation for repetition n
void CSimulator::outputSimulation(int n)
{

}



