/*
 * CSimulator.h
 *
 *  Created on: 26 Nov 2014
 *      Author: Jie Yang
 */

#ifndef CSIMULATOR_H_
#define CSIMULATOR_H_

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
	int nRepetitions, nYears, nOutputsPerYear, nTimeSteps;
	int nAG, CAGInfant, CAGPreSAC, CAGSAC, CAGAdult;
	int R0, lambda, LDecayRate;
	int TAGInfant, TAGPreSAC, TAGSAC, TAGAdult, treatStart, treatEnd, treatFreq;
	double dt;
	float gamma, z, k, sigma;
	float demog_eta, demog_b, InfantBeta, PreSACBeta, SACBeta, AdultBeta;
	float drugEfficacy, InfantCoverage, PreSACCoverage, SACCoverage, AdultCoverage;

	// General initialisation
	char buffer[BUFFER_SIZE];
	bool initialiseSimulation();
	CParamReader myReader;
	int nHosts;
	wormBurden** results; // Array of worm burdens
	CHost** hostPopulation; // Array of host population

	// Simulation
	void runSimulation();

	// Outputs
	void outputSimulation(int n);
};

#endif /* CSIMULATOR_H_ */
