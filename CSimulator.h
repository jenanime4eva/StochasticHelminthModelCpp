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
#include <random>


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
	double dt;


	// demographics.
	double demog_eta, demog_b, survivalDt;
	int survivalMaxIndex;
	double* survivalCurve;
	double* survivalCurveIntegral;
	double drawLifespan();

	/*
	int nAG, CAGInfant, CAGPreSAC, CAGSAC, CAGAdult;
	int R0, lambda, LDecayRate;
	int TAGInfant, TAGPreSAC, TAGSAC, TAGAdult, treatStart, treatEnd, treatFreq;
	float gamma, z, k, sigma;
	float InfantBeta, PreSACBeta, SACBeta, AdultBeta;
	float drugEfficacy, InfantCoverage, PreSACCoverage, SACCoverage, AdultCoverage;
	*/

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

	////////////////////////////////////////////////////
	/// Auxiliary functions.
	// a function to read in a list of doubles into a vector.
	double* readDoublesVector(char* valuesString, int& currentLength);
	// return an index for the value that is just above a uniform random deviate.
	int multiNomBasic(double* array, int length, double randNum);
	// a uniform random number generator (do better than this!)
	double myRand();
};

#endif /* CSIMULATOR_H_ */
