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

struct resultsCollection
{
	double time;
	double meanFemaleWorms;
	double meanInfantFemaleWorms;
	double meanPreSACFemaleWorms;
	double meanSACFemaleWorms;
	double meanAdultFemaleWorms;
};

class CSimulator {
public:
	CSimulator();
	virtual ~CSimulator();

	// File input/output
	ofstream logStream; // Output to log file
	bool initialiseIO(char* run, char* path, char* paramFilePath);
	int nRepetitions;
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
	double lambda, ReservoirDecayRate;

	// Treatment
	double* treatmentBreaks;
	int treatmentBreaksLength;
	double* coverage;
	int coverageLength;
	double drugEff;
	int treatStart;
	int treatEnd;
	double treatInterval;
	double* treatmentTimes;
	int treatmentTimesLength;

	// General initialisation
	char buffer[BUFFER_SIZE];
	bool initialiseSimulation();
	CParamReader myReader;
	int nHosts;
	double nYears;
	CRealization** myRealization; // Array of realisations

	// Simulation
	void runSimulation();

	// Outputs
	void outputSimulation();
	double surveyTimesDt;
	double* surveyResultTimes;
	int surveyResultTimesLength;
	double meanFemaleWorms;
	double meanInfantFemaleWorms;
	double meanPreSACFemaleWorms;
	double meanSACFemaleWorms;
	double meanAdultFemaleWorms;
	std::string resultsStub;

	// AUXILIARY FUNCTIONS

	// A function to read in a list of doubles into a vector
	double* readDoublesVector(char* valuesString, int& currentVectorLength);
	double* vectorArray;
	// Return an index for the value that is just above a uniform random deviate (for drawing lifespans)
	int multiNomBasic1(double* array, int length, double randNum);
	// Return an index for the value that is just above a uniform random deviate (for enacting an event in CRealization)
	int multiNomBasic2(double* array, int length, double randNum);
	// Calculate the psi value
	double calculatePsi();
	// Draw a lifespan from the population survival curve
	double drawLifespan();
	// Contact index array
	int* contactAgeGroupIndex();
	int* contactIndices;
	//Treatment index array
	int* treatmentAgeGroupIndex();
	int* treatmentIndices;
	// A uniform random number generator
	double myRandUniform();
	// A gamma distribution random number generator
	double myRandGamma(double l, double s);
	// A poisson distribution random number generator
	double myRandPoisson(double mu);
	// A binomial distribution random number generator
	double myRandBinomial(long n, double p);
	// An exponential distribution random number generator
	double myRandExponential(double a);
	// Find the sum of an array
	double sumArray(double* array,int arrayLength);
	// Find minimum of a list of values
	double min(double* Numbers, int Count);
	// Find maximum of a list of values
	double max(double* Numbers, int Count);
	// Return index of smallest element in an array
	int indexSmallestElement(double* array, int size);
};

#endif /* CSIMULATOR_H_ */
