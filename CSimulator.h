/*
 * File Name: CSimulator.h
 *
 * Created on: 26 Nov 2014
 * Author: Jie Yang
 *
 */

#ifndef CSIMULATOR_H_
#define CSIMULATOR_H_

#include <fstream>
#include <random>
#include <string>
#include "CParamReader.h"
#include "CHost.h"
#include "CRealization.h"

#define BUFFER_SIZE 1024 // Define buffer size

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
	bool initialiseIO(char* run, char* path, char* paramFilePath);
	int nRepetitions, nOutputsPerYear, nTimeSteps;
	double dt;
	std::string runName;
	std::string thePath;

	// Demographic structure
	double demogDt;
	double drawLifespan();
	double* hostMuData;  // Death rates
	int hostMuDataLength;
	double* muDataUpperBounds; // Upper bounds for death rates.
	int muDataUpperBoundsLength;
	int maxDtIntervals;
	double* survivalCurve;  // Survival to end of ith dt
	double* survivalCurveCumul;  // Cumulative sum of the above
	double* hostMu;
	double* probDeath; // Probability of dying in the ith dt
	double* probDeathIntegral; // Integral to end of ith dt
	double upperAgeBound;

	// Social structure
	double* contactAgeBreaks;
	int contactAgeBreaksLength;
	double* betaValues;
	int betaValuesLength;
	double* rhoValues;
	int rhoValuesLength;

	// Epidemiological parameters
	double k, R0, sigma, gamma;
	int lambda, ReservoirDecayRate;

	// Treatment
	double* treatmentBreaks;
	int treatmentBreaksLength;
	double* coverage;
	int coverageLength;
	double drugEff, treatInterval;
	int treatStart, nRounds;

	// General initialisation
	char buffer[BUFFER_SIZE];
	bool initialiseSimulation();
	CParamReader myReader;
	int nHosts;
	//wormBurden** results; // Array of worm burdens
	double nYears;
	CRealization myRealization; // DEBUG: just one at the moment

	// Simulation
	void runSimulation();

	// Outputs
	void outputSimulation();
	double* surveyResultTimes;
	int surveyResultTimesLength;
	std::string resultsStub;

	// AUXILIARY FUNCTIONS

	// A function to read in a list of doubles into a vector
	double* readDoublesVector(char* valuesString, int& currentLength);
	// Return an index for the value that is just above a uniform random deviate
	int multiNomBasic(double* array, int length, double randNum);
	// A uniform random number generator (do better than this!)
	double myRandUni();
};

#endif /* CSIMULATOR_H_ */
