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

using namespace std;

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
	ofstream logStream; // Output to log file
	bool initialiseIO(char* run, char* path, char* paramFilePath);
	int nRepetitions, nOutputsPerYear, nTimeSteps;
	double dt;
	string runName;
	string thePath;

	// Social structure
	double* contactAgeBreaks;
	int contactAgeBreaksLength;
	double* betaValues;
	int betaValuesLength;
	double* rhoValues;
	int rhoValuesLength;

	// Demographic structure
	double demogDt;
	double* hostMuData;  // Death rates
	int hostMuDataLength;
	double* muDataUpperBounds; // Upper bounds for death rates.
	int muDataUpperBoundsLength;
	int maxDtIntervals;
	int maxHostAge;
	double* maxHostAgeCompare;
	double* survivalCurve;  // Survival to end of ith dt
	double* survivalCurveCumul;  // Cumulative sum of the above
	double* hostMu;
	double* probDeath; // Probability of dying in the ith dt
	double* probDeathIntegral; // Integral to end of ith dt
	double upperAgeBound;
	double tinyIncrement;

	// Epidemiological parameters
	double k, R0, sigma, gamma, z, psi;
	int lambda, ReservoirDecayRate;

	// Treatment
	double* treatmentBreaks;
	int treatmentBreaksLength;
	double* coverage;
	int coverageLength;
	double drugEff, treatInterval;
	int treatStart, treatEnd;

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
	int surveyTimesBegin;
	int surveyTimesEnd;
	double surveyTimesDt;
	double* surveyResultTimes;
	int surveyResultTimesLength;
	std::string resultsStub;

	// AUXILIARY FUNCTIONS

	// A function to read in a list of doubles into a vector
	double* readDoublesVector(char* valuesString, int& currentVectorLength);
	double* vectorArray;
	// Return an index for the value that is just above a uniform random deviate
	int multiNomBasic(double* array, int length, double randNum);
	// Calculate the psi value
	double calculatePsi();
	// Draw a lifespan from the population survival curve
	double drawLifespan();
	// A uniform random number generator
	double myRandUniform();
	// A gamma distribution random number generator
	double myRandGamma();
	// A poisson distribution random number generator
	double myRandPoisson();
	// A binomial distribution random number generator
	double myRandBinomial();
	// Find the cumsum of a vector
	//vector<double> cumsum(vector<double> x);
	// Find the sum of a vector
	//double sum(vector<double> x);
	// Find minimum of a list of values
	double min(double* Numbers, int Count);
	// Find maximum of a list of values
	double max(double* Numbers, int Count);
};

#endif /* CSIMULATOR_H_ */
