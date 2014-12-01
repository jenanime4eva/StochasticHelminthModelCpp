/*
 * CSimulator.h
 *
 *  Created on: 26 Nov 2014
 *      Author: Jie Yang
 */

#pragma once

#include <fstream>
#include "CParamReader.h"
#include "CHost.h"

#define BUFFER_SIZE 1024

struct wormBurden
{
	double time;
	int nWorms;
};

class CSimulator {
public:
	CSimulator();
	virtual ~CSimulator();

	// File input/output
	std::ofstream logStream, resultsStream; // Output to log file and results file
	bool initialiseIO(char* logFileName, char* paramFileName, char* resultsFileName);

	// Initialisation
	char buffer[BUFFER_SIZE];
	bool initialiseSimulation();
	CParamReader myReader;

	// Simulation
	int nRepetitions, nYears, nOutputsPerYear, nHosts;
	wormBurden** results; // Array of worm burdens
	CHost** hostPopulation;
	void runSimulation();

	// Outputs
	void outputSimulation(int n);
};
