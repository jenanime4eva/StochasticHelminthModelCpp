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
#include <vector.h>
#include "randlib.h"

using namespace std;

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

	// Social structure
	contactAgeBreaksLength = 0;
	betaValuesLength = 0;
	rhoValuesLength = 0;

	// Demographic structure
	survivalCurve = survivalCurveCumul = NULL;
	hostMu = probDeath = probDeathIntegral = NULL;
	demogDt = 0;
	hostMuDataLength = 0;
	muDataUpperBoundsLength = 0;
	upperAgeBound = 0;
	maxDtIntervals = 0;

	// Epidemiological parameters
	k = 0;
	lambda = 0;
	R0 = 0;
	ReservoirDecayRate = 0;
	sigma = 0;
	gamma = 0;

	// Treatment parameters
	treatmentBreaksLength = 0;
	coverageLength = 0;
	drugEff = 0;
	treatStart = 0;
	nRounds = 0;
	treatInterval = 0;

	// Results
	surveyTimesBegin = 0;
	surveyTimesEnd = 0;
	surveyTimesDt = 0;
	surveyResultTimes = NULL;
	surveyResultTimesLength = 0;
}

// Class Destructor
CSimulator::~CSimulator() {
	if (survivalCurve != NULL)
		delete[] survivalCurve;

	if (survivalCurveCumul != NULL)
		delete[] survivalCurveCumul;

	if (hostMu != NULL)
		delete[] hostMu;

	if (probDeath != NULL)
		delete[] probDeath;

	if (probDeathIntegral != NULL)
		delete[] probDeathIntegral;

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

	logStream << "\nTEST PARAMETERS READ IN\n\n" << std::flush;

	// Number of repetitions
	nRepetitions = atoi(myReader.getParamString("repNum"));
	logStream << "Number of repetitions: " << nRepetitions << "\n" << std::flush; // Test flag

	// Number of years to run
	nYears = atoi(myReader.getParamString("nYears"));
	logStream << "Number of years to run: " << nYears << "\n" << std::flush; // Test flag

	// Number of hosts
	nHosts = atoi(myReader.getParamString("nHosts"));
	logStream << "Number of hosts: " << nHosts << "\n" << std::flush; // Test flag

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

	/// SET UP DEMOGRAPHY

	// Read in host death rates
	temp = myReader.getParamString("hostMuData");
	if (temp != NULL) {
		hostMuData = readDoublesVector(temp, hostMuDataLength);
	}

	// Read in host death rate upper bounds
	temp = myReader.getParamString("upperBoundData");
	if (temp != NULL) {
		muDataUpperBounds = readDoublesVector(temp, muDataUpperBoundsLength);
	}

	// Age step for survival curve in years
	demogDt = atof(myReader.getParamString("demogDt"));

	// Construct mu, probability of death and survival vectors
	int currentMuIndex = 0;
	int maxHostAge = 0;
	double currentSurvival = 1;
	double currentSurvivalCumul = 0;
	double currentProbDeathCumul = 0;
	double currentMuDtCumul = 0;
	double tinyIncrement = 0.01;

	maxDtIntervals = (int) floor(muDataUpperBounds[muDataUpperBoundsLength-1]/demogDt);
	upperAgeBound = maxDtIntervals * demogDt;

	vector<double> maxHostAgeCompare = {upperAgeBound,contactAgeBreaks[contactAgeBreaksLength-1]}; // Create comparison array
	int maxHostAgeCompareLength = maxHostAgeCompare.size(); // Get length of list of values
	maxHostAge = (int) min(maxHostAgeCompare,maxHostAgeCompareLength); 	// Get maximum host age

	// Recalculate maxDtIntervals now that we know the maxHostAge value
	maxDtIntervals = maxHostAge/demogDt;

	survivalCurve = new double[maxDtIntervals];
	survivalCurveCumul = new double[maxDtIntervals];
	hostMu = new double[maxDtIntervals];
	probDeath = new double[maxDtIntervals];
	probDeathIntegral = new double[maxDtIntervals];

	for (int i=0;i<maxDtIntervals;i++)
	{
		double currentIntEnd = (i + 1) * demogDt;
		// Is current dt interval within data upper bound?
		if (muDataUpperBounds[currentMuIndex] + tinyIncrement < currentIntEnd)
		{
			currentMuIndex++; // Add one to currentMuIndex
		}
		hostMu[i] = hostMuData[currentMuIndex];

		probDeath[i] = currentSurvival * hostMu[i] * demogDt;
		probDeathIntegral[i] = currentProbDeathCumul + probDeath[i];
		currentProbDeathCumul = probDeathIntegral[i]; // Cumulative probability of dying in the ith year

		currentMuDtCumul += hostMu[i] * demogDt;
		survivalCurve[i] = exp(-currentMuDtCumul);
		currentSurvival = survivalCurve[i]; // Host survival curve

		survivalCurveCumul[i] = survivalCurve[i] + currentSurvivalCumul;
		currentSurvivalCumul = survivalCurveCumul[i];
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
	double z = exp(-gamma); // Fecundity parameter

	// Dummy psi value prior to R0 calculation
	double psi = 1.0;

	// SET UP TREATMENT

	// Read in treatment age group breaks
	temp = myReader.getParamString("treatmentBreaks");
	if (temp != NULL) {
		treatmentBreaks = readDoublesVector(temp, treatmentBreaksLength);
	}
	logStream << "\ntreatmentBreaks vector length: " << treatmentBreaksLength << "\n" << std::flush; // Test flag
	logStream << "Infant treatment break: " << treatmentBreaks[0] << "\n" << std::flush; // Test flag
	logStream << "pre-SAC treatment break: " << treatmentBreaks[1] << "\n" << std::flush; // Test flag
	logStream << "SAC treatment break: " << treatmentBreaks[2] << "\n" << std::flush; // Test flag
	logStream << "Adult treatment break: " << treatmentBreaks[3] << "\n" << std::flush; // Test flag
	logStream << "Maximum treatment age: " << treatmentBreaks[4] << "\n" << std::flush; // Test flag

	// Read in coverages
	temp = myReader.getParamString("coverage");
	if (temp != NULL) {
		coverage = readDoublesVector(temp, coverageLength);
	}
	logStream << "\nCoverages vector length: " << coverageLength << "\n" << std::flush; // Test flag
	logStream << "Infant coverage: " << coverage[0] << "\n" << std::flush; // Test flag
	logStream << "pre-SAC coverage: " << coverage[1] << "\n" << std::flush; // Test flag
	logStream << "SAC coverage: " << coverage[2] << "\n" << std::flush; // Test flag
	logStream << "Adult coverage: " << coverage[3] << "\n\n" << std::flush; // Test flag

	// Drug efficacy
	drugEff = atof(myReader.getParamString("drugEff"));

	// Treatment year start
	treatStart = atoi(myReader.getParamString("treatStart"));

	// Number of treatment rounds
	nRounds = atoi(myReader.getParamString("nRounds"));

	// Interval between treatments in years
	treatInterval = atof(myReader.getParamString("treatInterval"));

	// SET UP RESULTS COLLECTION

	// Year to start getting data from all individuals in the population
	surveyTimesBegin = atoi(myReader.getParamString("surveyTimesBegin"));
	// Year to stop getting data from all individuals in the population
	surveyTimesEnd = atoi(myReader.getParamString("surveyTimesEnd"));
	// Time step for the survey times in years
	surveyTimesDt = atof(myReader.getParamString("surveyTimesDt"));
	// Now create a vector of survey result times
	surveyResultTimesLength = ((surveyTimesEnd - surveyTimesBegin)/surveyTimesDt)+1;
	surveyResultTimes = new double[surveyResultTimesLength];
	surveyResultTimes[0] = surveyTimesBegin; // First vector entry
	for (int i=1;i<surveyResultTimesLength;i++)
	{
		surveyResultTimes[i] = surveyResultTimes[i-1] + surveyTimesDt;
	}
	logStream << "surveyResultTimesLength: " << surveyResultTimesLength << "\n" << std::flush; // Test flag
	logStream << "surveyResultTimes first entry: " << surveyResultTimes[0] << "\n" << std::flush; // Test flag
	logStream << "surveyResultTimes second entry: " << surveyResultTimes[1] << "\n" << std::flush; // Test flag
	logStream << "surveyResultTimes last entry: " << surveyResultTimes[surveyResultTimesLength-1] << "\n" << std::flush; // Test flag


	// Set up a realisation
	myRealization.initialize(this);
	myRealization.run();

	/*// Leave this commented for now
	// Delete memory related to the variables using the readDoublesVector function to avoid memory leaks
	if (contactAgeBreaks != NULL)
		delete[] contactAgeBreaks;

	if (betaValues != NULL)
		delete[] betaValues;

	if (rhoValues != NULL)
			delete[] rhoValues;

	if (hostMuData != NULL)
		delete[] hostMuData;

	if (muDataUpperBounds != NULL)
		delete[] muDataUpperBounds;

	if (treatmentBreaks != NULL)
		delete[] treatmentBreaks;

	if (coverage != NULL)
		delete[] coverage;
	*/

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

// Creates a vector of doubles from string from the parameter file
vector<double> CSimulator::readDoublesVector(char* currentString, int& currentVectorLength)
{
	// IMPORTANT: Make sure the string read in does not have trailing zeros
	// These are currently removed in the parameter reading function

	char* endPointer; // endPointer for each call to strto_ type functions

	// This part is just to count the number of entries in the vector
	int counter = 1;
	double* tempVector = new double[strlen(currentString)]; // Create temporary vector
	tempVector[0] = strtod(currentString, &endPointer);
	while(strlen(endPointer)>0)
	{
		tempVector[counter] = strtod(endPointer, &endPointer);
		counter++; // Add one to counter
	}
	// Delete tempVector as we don't need it anymore
	if (tempVector!=NULL)
		delete[] tempVector;

	// NOW create the vector of doubles of values we actually require
	vector<double> vectorArray(counter);
	vectorArray[0] = strtod(currentString, &endPointer);
	for(int i=1;i<counter;i++)
	{
		vectorArray[i] = strtod(endPointer, &endPointer);
	}

	currentVectorLength = vectorArray.size(); // Count number of elements in the vector

	return vectorArray;
}

// This function takes a random number (0-1) and multiplies it by the max of the passed array.
// It then finds the smallest index that has an array value greater than the product above.
// For a cumulative multinomial array, this will return the index of the event that occurred.
// Also used for drawing a lifespan from the survival curve integral.
int CSimulator::multiNomBasic(double* array, int length, double randNum)
{
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

// Calculate the psi value
double CSimulator::calculatePsi()
{
	double psi = 1;
	return 2*psi; // This 2 is here because it's infection with both male and female worms and only half of these are going to be female
}

// Draw a life span from the survival curve from the population
double CSimulator::drawLifespan()
{
	// Get a random integer from the probDeathIntegral using the multinomial generator. This shouldn't be zero!
	double currentRand = myRandUni();

	int index = multiNomBasic(probDeathIntegral, maxDtIntervals, currentRand); // maxDtIntervals = maxHostAge here?

	// Choose a point in the middle of the interval in which the person dies
	double ans = (index + 0.5) * demogDt;

	return ans;
}

// Uniform random number generator (should be using the <random> functions)
double CSimulator::myRandUni()
{
	//genunf(double low, double high) generates uniform real between low and high
	  return  genunf(0,1);

}

// Function to find minimum of a list of values
double CSimulator::min(vector<double> Numbers, int Count)
{
	double Minimum = Numbers[0];

	for(int i=0;i<Count;i++)
	{
		if(Minimum>Numbers[i])
		{
			Minimum = Numbers[i];
		}
	}

	return Minimum;
}

// Function to find maximum of a list of values
double CSimulator::max(vector<double> Numbers, int Count)
{
	double Maximum = Numbers[0];

	for(int i=0;i<Count;i++)
	{
		if(Maximum<Numbers[i])
		{
			Maximum = Numbers[i];
		}
	}

	return Maximum;
}
