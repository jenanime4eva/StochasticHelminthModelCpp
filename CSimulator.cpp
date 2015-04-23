/*
 * File Name: CSimulator.cpp
 *
 *  Created on: 26 Nov 2014
 *  Author: Jie Yang
 *
 */

#include "CSimulator.h"
#include "CHost.h"
#include "CParamReader.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "randlib.h"

// Class Constructor
CSimulator::CSimulator() {
	// Set default values

	// Model running parameters
	nRepetitions = 0;
	nYears = 0;
	nHosts = 0;
	nOutputsPerYear = 0;
	nTimeSteps = 0;
	dt = 0;

	// Demographic structure
	survivalCurve = survivalCurveCumul = NULL;
	hostMu = probDeath = probDeathIntegral = NULL;
	demogDt = 0;
	hostMuData = NULL;
	hostMuDataLength = 0;
	muDataUpperBounds = NULL;
	muDataUpperBoundsLength = 0;
	upperAgeBound = 0;
	maxDtIntervals = 0;

	// Social structure
	contactAgeBreaks = NULL;
	contactAgeBreaksLength = 0;
	betaValues = NULL;
	betaValuesLength = 0;
	rhoValues = NULL;
	rhoValuesLength = 0;

	// Epidemiological parameters
	k = 0;
	lambda = 0;
	R0 = 0;
	ReservoirDecayRate = 0;
	sigma = 0;
	gamma = 0;

	// Treatment parameters
	treatmentBreaks = NULL;
	treatmentBreaksLength = 0;
	coverage = NULL;
	coverageLength = 0;
	drugEff = 0;
	treatStart = 0;
	nRounds = 0;
	treatInterval = 0;

	// Results
	surveyResultTimes = NULL;
	surveyResultTimesLength = 0;
}

// Class Destructor
CSimulator::~CSimulator() {
	if (survivalCurve != NULL)
		delete[] survivalCurve;

	if (survivalCurveCumul != NULL)
		delete[] survivalCurveCumul;

	if (hostMuData != NULL)
		delete[] hostMuData;

	if (muDataUpperBounds != NULL)
		delete[] muDataUpperBounds;

	if (hostMu != NULL)
		delete[] hostMu;

	if (probDeath != NULL)
		delete[] probDeath;

	if (probDeathIntegral != NULL)
		delete[] probDeathIntegral;

	if (contactAgeBreaks != NULL)
		delete[] contactAgeBreaks;

	if (betaValues != NULL)
		delete[] betaValues;

	if (rhoValues != NULL)
		delete[] rhoValues;

	if (treatmentBreaks != NULL)
		delete[] treatmentBreaks;

	if (coverage != NULL)
		delete[] coverage;

	if (surveyResultTimes != NULL)
		delete[] surveyResultTimes;
}

// Initialise the input/output aspects of the simulator
bool CSimulator::initialiseIO(char* run, char* path, char* paramFilePath)
{
	// As string classes
	runName = run;
	thePath = path;

	// Setting up all the files needed
	std::string logFilePath = thePath + runName + ".log.txt";
	logStream.open(logFilePath.c_str());
	if (!logStream.is_open()) {
		std::cout << "Couldn't open the log file: " << logFilePath
				<< "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Log file opened.\n" << std::flush;

	myReader.setNewFileName(paramFilePath);
	if (!myReader.setNewFileName(paramFilePath)) {
		logStream << "Couldn't open the parameter file: " << paramFilePath
				<< "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Param reader configured.\n" << std::flush;

	// Stub name for all results files
	resultsStub = thePath + runName;
	logStream << "Stub name for all results files: " << resultsStub << "\n" << std::flush;

	// Everything is ok
	return true;
}

// General initialisation
bool CSimulator::initialiseSimulation()
{
	char* temp; // General purpose pointer to string

	// READ IN MODEL RUNNING PARAMETERS

	logStream << "\nPARAMETERS READ IN\n" << std::flush;

	// Number of repetitions
	nRepetitions = atoi(myReader.getParamString("repNum"));
	logStream << "Number of repetitions: " << nRepetitions << "\n" << std::flush; // Test flag

	// Number of years to run
	nYears = atoi(myReader.getParamString("nYears"));
	logStream << "Number of years to run: " << nYears << "\n" << std::flush; // Test flag

	// Number of hosts
	nHosts = atoi(myReader.getParamString("nHosts"));
	logStream << "Number of hosts: " << nHosts << "\n" << std::flush; // Test flag

	// SET UP DEMOGRAPHY

	// Read in host death rates
	temp = myReader.getParamString("hostMu");
	if (temp != NULL) {
		hostMuData = readDoublesVector(temp, hostMuDataLength);
	}
	logStream << "hostMuData vector length: " << hostMuDataLength << "\n" << std::flush; // Test flag

	// Read in host death rate upper bounds
	temp = myReader.getParamString("upperBoundData");
	if (temp != NULL) {
		muDataUpperBounds = readDoublesVector(temp, muDataUpperBoundsLength);
	}

	// Age step for survival curve in years
	demogDt = atof(myReader.getParamString("demogDt"));

	// Construct mu, probability of death and survival vectors
	maxDtIntervals = (int) floor(
			muDataUpperBounds[muDataUpperBoundsLength - 1] / demogDt);
	upperAgeBound = maxDtIntervals * demogDt;
	int currentMuIndex = 0;
	double currentSurvival = 1;
	double currentSurvivalCumul = 0;
	double currentProbDeathCumul = 0;
	double currentMuDtCumul = 0;
	double tinyIncrement = 0.01;
	survivalCurve = new double[maxDtIntervals];
	survivalCurveCumul = new double[maxDtIntervals];
	hostMu = new double[maxDtIntervals];
	probDeath = new double[maxDtIntervals];
	probDeathIntegral = new double[maxDtIntervals];

	for (int i = 0; i < maxDtIntervals; i++) {
		double currentIntEnd = (i + 1) * demogDt;
		// Is current dt interval within data upper bound?
		if (muDataUpperBounds[currentMuIndex] + tinyIncrement < currentIntEnd)
			currentMuIndex++;
		hostMu[i] = hostMuData[currentMuIndex];
		probDeath[i] = currentSurvival * hostMu[i] * demogDt;
		currentMuDtCumul += hostMu[i] * demogDt;

		probDeathIntegral[i] = currentProbDeathCumul + probDeath[i];
		currentProbDeathCumul = probDeathIntegral[i];

		survivalCurve[i] = exp(-currentMuDtCumul);
		currentSurvival = survivalCurve[i];

		survivalCurveCumul[i] = survivalCurve[i] + currentSurvivalCumul;
		currentSurvivalCumul = survivalCurveCumul[i];
	}

	// SET UP SOCIAL STRUCTURE

	// Read in contact age group breaks
	temp = myReader.getParamString("contactAgeBreaks");
	if (temp != NULL) {
		contactAgeBreaks = readDoublesVector(temp, contactAgeBreaksLength);
	}

	// Read in beta values (contact rates)
	temp = myReader.getParamString("betaValues");
	if (temp != NULL) {
		betaValues = readDoublesVector(temp, betaValuesLength);
	}

	// Read in rho values (contribution to the reservoir by contact age group)
	temp = myReader.getParamString("rhoValues");
	if (temp != NULL) {
		rhoValues = readDoublesVector(temp, rhoValuesLength);
	}

	// READ IN EPIDEMIOLOGICAL PARAMETERS

	// Shape parameter of assumed negative binomial distribution of worms amongst host
	k = atof(myReader.getParamString("k"));

	// Eggs per gram
	lambda = atoi(myReader.getParamString("lambda"));

	// Basic reproductive number
	R0 = atof(myReader.getParamString("R0"));
	logStream << "R0: " << R0 << "\n" << std::flush; // Test flag

	// Decay rate of eggs in the environment
	ReservoirDecayRate = atoi(myReader.getParamString("ReservoirDecayRate"));

	// Worm death rate i.e. 1/worm life span, same for all development stages
	sigma = atof(myReader.getParamString("sigma"));

	// Exponential density dependence of parasite adult stage (N.B. fecundity parameter z = exp(-gamma))
	gamma = atof(myReader.getParamString("gamma"));

	// SET UP TREATMENT

	// Read in treatment age group breaks
	temp = myReader.getParamString("treatmentBreaks");
	if (temp != NULL) {
		treatmentBreaks = readDoublesVector(temp, treatmentBreaksLength);
	}

	// Read in coverages
	temp = myReader.getParamString("coverage");
	if (temp != NULL) {
		coverage = readDoublesVector(temp, coverageLength);
	}

	// Drug efficacy
	drugEff = atof(myReader.getParamString("drugEff"));

	// Treatment year start
	treatStart = atoi(myReader.getParamString("treatStart"));

	// Number of treatment rounds
	nRounds = atoi(myReader.getParamString("nRounds"));

	// Interval between treatments in years
	treatInterval = atof(myReader.getParamString("treatInterval"));

	// SET UP RESULTS COLLECTION
	temp = myReader.getParamString("surveyTimes");
	if (temp != NULL) {
		surveyResultTimes = readDoublesVector(temp, surveyResultTimesLength);
	}

	// Set up a realisation
	myRealization.initialize(this);
	myRealization.run();

	return true;
}

// Run simulation
void CSimulator::runSimulation()
{

}

// What simulation outputs do we want to see?
// 1) Frequency vs Number of worms
// 2) Individual runs of worm burden across time (FOCUS ON THIS ONE FIRST)
// 3) Mean worm burden across time
// 4) Prevalence across time

// Output simulation
void CSimulator::outputSimulation()
{
	// Output survey-type results if there's any
	if (surveyResultTimes != NULL) {
		// There's some survey results
		// Create output stream
		std::string surveyResultsOut = thePath + runName + ".surveyResults.txt";
		std::ofstream surveyStream(surveyResultsOut.c_str());

		// DEBUG: currently only ONE realization, so...
		for (int j = 0; j < surveyResultTimesLength; j++)
			surveyStream << surveyResultTimes[j] << "\t";
		surveyStream << "\n";
		// Loop through the hosts
		for (int i = 0; i < nHosts; i++) {
			for (int j = 0; j < surveyResultTimesLength; j++) {
				surveyStream << myRealization.surveyResultsArray[j][i].age
						<< "\t";
			}
			surveyStream << "\n";
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
/// Auxiliary function definitions

// Creates a vector of doubles from string from param file
// Allocates memory so need to delete
double* CSimulator::readDoublesVector(char* currentString, int& currentLength)
{
	char* endPointer; // endPointer for each call to strto_ type functions

	// Read in the length of the vector
	int length = strtol(currentString, &endPointer, 10);
	currentLength = length;

	if (length < 0)
		return NULL;

	// Create array of doubles
	double* vectorArray = new double[length];
	for (int i = 0; i < length; i++) {
		vectorArray[i] = strtod(endPointer, &endPointer);
	}

	return vectorArray;
}

// This function takes a random number (0-1) and multiplies it by the max of the passed array.
// It then finds the smallest index that has an array value greater than the product above.
// For a cumulative multinomial array, this will return the index of the event that occurred.
// Also used for drawing a lifespan from the survival curve integral.
int CSimulator::multiNomBasic(double* array, int length, double randNum) {
	int loopMax = ceil(log(length) / log(2) + 2);
	int bottom = -1;
	//double bottomVal;
	int top = length - 1;
	double topVal = array[top];
	double target = topVal * randNum;
	int count = 0;

	while (++count < loopMax && (top - bottom > 1)) {
		int mid = (top + bottom) / 2;
		double midVal = array[mid];
		if (midVal >= target) {
			top = mid;
			topVal = midVal;
		} else {
			bottom = mid;
			//bottomVal = midVal;
		}
	}

	if (count >= loopMax) {
		logStream << "Max iterations exceeded in multiNomBasic(...),\n"
				<< std::flush;
		return -1;
	}

	return top;
}

// Uniform random number generator (should be using the <random> functions)
double CSimulator::myRandUni()
{
	//genunf(double low, double high) generates uniform real between low and high
	  return  genunf(0,1);

}

// Draw a life span from the survival curve from the population
double CSimulator::drawLifespan() {
	// Get a random integer from the probDeathIntegral using the multinomial generator. This shouldn't be zero!
	double currentRand = myRandUni();

	int index = multiNomBasic(probDeathIntegral, maxDtIntervals, currentRand);

	// Choose a point in the middle of the interval in which the person dies
	double ans = (index + 0.5) * demogDt;

	return ans;
}
