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
	vector<double> contactAgeBreaks;
	int contactAgeBreaksLength;
	vector<double> betaValues;
	int betaValuesLength;
	vector<double> rhoValues;
	int rhoValuesLength;

	// Demographic structure
	double demogDt;
	double drawLifespan();
	vector<double> hostMuData;  // Death rates
	int hostMuDataLength;
	vector<double> muDataUpperBounds; // Upper bounds for death rates.
	int muDataUpperBoundsLength;
	int maxDtIntervals;
	double* survivalCurve;  // Survival to end of ith dt
	double* survivalCurveCumul;  // Cumulative sum of the above
	double* hostMu;
	double* probDeath; // Probability of dying in the ith dt
	double* probDeathIntegral; // Integral to end of ith dt
	double upperAgeBound;

	// Epidemiological parameters
	double k, R0, sigma, gamma;
	int lambda, ReservoirDecayRate;

	// Treatment
	vector<double> treatmentBreaks;
	int treatmentBreaksLength;
	vector<double> coverage;
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
	int surveyTimesBegin;
	int surveyTimesEnd;
	double surveyTimesDt;
	double* surveyResultTimes;
	int surveyResultTimesLength;
	std::string resultsStub;

	// AUXILIARY FUNCTIONS

	// A function to read in a list of doubles into a vector
	vector<double> readDoublesVector(char* valuesString, int& currentVectorLength);
	// Return an index for the value that is just above a uniform random deviate
	int multiNomBasic(double* array, int length, double randNum);
	// A uniform random number generator (do better than this!)
	double myRandUni();
	// Find minimum of a list of values
	double min(double* Numbers, int Count);
	// Find maximum of a list of values
	double max(double* Numbers, int Count);
};

#endif /* CSIMULATOR_H_ */
