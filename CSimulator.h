/*
 * CSimulator.h
 *
 *  Created on: 26 Nov 2014
 *      Author: Jie Yang
 */

#ifndef CSIMULATOR_H_
#define CSIMULATOR_H_

#include <fstream>
#include <random>
#include <string>

#include "CParamReader.h"
#include "CHost.h"
#include "CRealization.h"

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
	std::ofstream logStream; // Output to log file
	bool initialiseIO(char* logFileName, char* paramFileName, char* resultsFileName);
	int nRepetitions, nOutputsPerYear, nTimeSteps;
	double dt;
	std::string runName;
	std::string thePath;


	// demographics.
	double demog_eta, demog_b;
	double demogDt;
	double drawLifespan();
	double* hostMuData;  // death rates.
	int hostMuDataLength;
	double* muDataUpperBounds; // upper bounds for death rates.
	int muUpperBoundsLength;

	int maxDtIntevals;
	double* survivalCurve;  // survival to end of ith dt.
	double* survivalCurveCumul;  // cumulative sum of the above.
	double* hostMu;
	double* probDeath; // prob of dying in the ith dt.
	double* probDeathIntegral; // integral to end of ith dt.
	double upperAgeBound;

	// Results.
	double* surveyResultTimes; // default NULL. If NULL none are collected.
	int surveyResultTimesLength;
	std::string resultsStub;

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
	//wormBurden** results; // Array of worm burdens
	double startYear, nYears;
	CRealization myRealization; // DEBUG: just one at the moment.

	// Simulation
	void runSimulation();

	// Outputs
	void outputSimulation();

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
