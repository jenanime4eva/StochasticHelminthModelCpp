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
//#include <random>

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


	// Demographics
	double demog_eta, demog_b, survivalDt,survivalMaxAge;
	int survivalMaxIndex;
	double* survivalCurve;
	double* survivalCurveIntegral;
	double drawLifespan();

	int CAGInfant, CAGPreSAC, CAGSAC, CAGAdult;
	int R0, lambda, LDecayRate;
	int TAGInfant, TAGPreSAC, TAGSAC, TAGAdult, treatStart, treatEnd, treatFreq;
	double gamma, z, k, sigma;
	double InfantBeta, PreSACBeta, SACBeta, AdultBeta;
	double drugEfficacy, InfantCoverage, PreSACCoverage, SACCoverage, AdultCoverage;

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

	////////////////////////////////////////////////////////////////////////////////
	/// Auxiliary functions
	char* endPointer;
	// A function to read in a list of doubles into a vector
	double* readDoublesVector(char* valuesString);
	// A function to read in a list of integers into a vector
	int* readIntsVector(char* valuesString);
	// Return an index for the value that is just above a uniform random deviate
	int multiNomBasic(double* array, int length, double randNum);
//	// A uniform random number generator (do better than this!) -> GET RANDOM HEADER WORKING HERE
//	double myRand();
};

#endif /* CSIMULATOR_H_ */
