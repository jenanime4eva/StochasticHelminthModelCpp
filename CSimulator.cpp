/*
 * CSimulator.cpp
 *
 *  Created on: 26 Nov 2014
 *      Author: Jie Yang
 */

#include ".\CSimulator.h"
#include "CHost.h"
#include "CParamReader.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "randlib.h"

// Class Constructor
CSimulator::CSimulator()
{
	// Default values
	hostPopulation = NULL;
	results = NULL;
	nRepetitions = 0;
	nYears = 0;
	nOutputsPerYear = 0;
	nTimeSteps = 0;
	dt = 0;

	nHosts = 0;
	demog_b = demog_eta = 0;
	survivalCurve = survivalCurveIntegral = NULL;
	survivalDt = 0;
	survivalMaxAge = 0;
	survivalMaxIndex = 0;

	CAGInfant = CAGPreSAC = CAGSAC = CAGAdult = 0;
	InfantBeta = PreSACBeta = SACBeta = AdultBeta = 0;

	R0 = 0;
	lambda = 0;
	gamma = 0;
	z = 0;
	k = 0;
	sigma = 0;
	LDecayRate = 0;

	TAGInfant = TAGPreSAC = TAGSAC = TAGAdult = 0;
	InfantCoverage = PreSACCoverage = SACCoverage = AdultCoverage = 0;
	drugEfficacy = 0;
	treatStart = treatEnd = treatFreq = 0;

	endPointer = NULL;
}

// Class Destructor
CSimulator::~CSimulator()
{
	// Delete results array allocations
	if (results!=NULL)
	{
		// Delete memory
		for (int i=0;i<nRepetitions;i++)
		{
			delete results[i];
		}
		delete[] results;
	}

	// Delete hostPopulation array allocations
	if (hostPopulation!=NULL)
	{
		for (int i=0;i<nHosts;i++)
		{
			delete hostPopulation[i];
		}
		delete[] hostPopulation;
	}

	if(survivalCurve!=NULL)
		delete[] survivalCurve;

	if(survivalCurveIntegral!=NULL)
		delete[] survivalCurveIntegral;

}

// Initialise the input/output aspects of the simulator
bool CSimulator::initialiseIO(char* logFileName, char* paramFileName, char* resultsFileName)
{
	// Setting up all the files needed
	logStream.open(logFileName);
	if(!logStream.is_open())
	{
		std::cout << "Couldn't open the log file: " << logFileName << "\nexiting\n";
		return false;
	}
	logStream << "Log file opened.\n" << std::flush;

	if (!myReader.setNewFileName(paramFileName))
 	{
		logStream << "Couldn't open the param file: " << paramFileName << "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Param reader configured.\n" << std::flush;

	resultsStream.open(resultsFileName);
	if(!resultsStream.is_open())
	{
		logStream << "Couldn't open the results file: " << resultsFileName << "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Results file opened.\n" << std::flush;

	// Everything is ok
	return true;
}

// General initialisation
bool CSimulator::initialiseSimulation()
{
	// Get model related parameters

	int vectorLength;

	// Number of repetitions
	nRepetitions = atoi(myReader.getParamString("numberRepetitions"));

	// Number of years to run
	nYears = atoi(myReader.getParamString("numberYears"));

	// Number of outputs per year
	nOutputsPerYear = atoi(myReader.getParamString("numberOutputsPerYear"));

	// Number of hosts
	nHosts = atoi(myReader.getParamString("numberHosts"));

	///////////////////////////////////////////////////////////////////////////////////////////////
	///  Set up demography

	// Read a parameter to determine whether to use data for the survival curve or a named function

	// Use the exponential-exponential function to define the survival curve
	double* demography = readDoublesVector(myReader.getParamString("demography"),vectorLength);
	demog_eta = demography[0];
	demog_b = demography[1];

	survivalDt = atof(myReader.getParamString("survivalDt"));
	survivalMaxAge = atof(myReader.getParamString("survivalMaxAge"));

	survivalMaxIndex = (int) ceil(survivalMaxAge/survivalDt);
	survivalCurve = new double[survivalMaxIndex];
	survivalCurveIntegral = new double[survivalMaxIndex];
	double subTotal = 0;
	for(int i=0;i<survivalMaxIndex;i++)
	{
		// Exponential-exponential function
		// survivalCurve[i] = exp(-demog_eta*(exp(demog_b*survivalDt*i)-1)); // Uncomment later
		survivalCurve[i] = exp(-demog_eta*survivalDt*i);
		subTotal += survivalCurve[i]; // subTotal = subTotal + survivalCurve
		survivalCurveIntegral[i] = (subTotal - (survivalCurve[0]+survivalCurve[i])/2)*survivalDt;
	}

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;
		// Do lifespan stuff...
		double lifespan = drawLifespan();
		hostPopulation[i] -> birthDate = -genunf(0,1)*lifespan;
		hostPopulation[i] -> deathDate = hostPopulation[i] -> birthDate + lifespan;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	//// DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE ///

	// Contact age groups
	int* contactAgeGroups = readIntsVector(myReader.getParamString("contactAgeGroups"),vectorLength);
	CAGInfant = contactAgeGroups[0];
	CAGPreSAC = contactAgeGroups[1];
	CAGSAC = contactAgeGroups[2];
	CAGAdult = contactAgeGroups[3];

	// Beta values
	double* beta = readDoublesVector(myReader.getParamString("beta"),vectorLength);
	InfantBeta = beta[0];
	PreSACBeta = beta[1];
	SACBeta = beta[2];
	AdultBeta = beta[3];

	// R0
	R0 = atoi(myReader.getParamString("R0"));

	// Egg per gram
	lambda = atoi(myReader.getParamString("lambda"));

	// Exponential density dependence of parasite adult stage
	gamma = atof(myReader.getParamString("gamma"));

	// Fecundity parameter z = exp(-gamma)
	z = atof(myReader.getParamString("z"));

	// Shape parameter of assumed negative binomial distribution of worms amongst host
	k = atof(myReader.getParamString("k"));

	// Worm death rate i.e. 1/worm_life_span
	sigma = atof(myReader.getParamString("sigma"));

	// Reservoir decay rate (decay rate of eggs in the environment)
	LDecayRate = atoi(myReader.getParamString("LDecayRate"));

	// Treatment age groups
	int* treatmentAgeGroups = readIntsVector(myReader.getParamString("treatmentAgeGroups"),vectorLength);
	TAGInfant = treatmentAgeGroups[0];
	TAGPreSAC = treatmentAgeGroups[1];
	TAGSAC = treatmentAgeGroups[2];
	TAGAdult = treatmentAgeGroups[3];

	// Coverages
	double* coverages = readDoublesVector(myReader.getParamString("coverages"),vectorLength);
	InfantCoverage = coverages[0];
	PreSACCoverage = coverages[1];
	SACCoverage = coverages[2];
	AdultCoverage = coverages[3];

	// Drug efficacy
	drugEfficacy = atof(myReader.getParamString("drugEfficacy"));

	// Chemotherapy timings
	int* treatmentTimes = readIntsVector(myReader.getParamString("treatmentTimes"),vectorLength);
	treatStart = treatmentTimes[0];
	treatEnd = treatmentTimes[1];
	treatFreq = treatmentTimes[2];

	// Print some parameters of different types to log file to check if they are being read in correctly
	logStream << "\nNumber of repetitions (int type): " << nRepetitions << "\n"
				<< "Fecundity parameter (double type): " << z << "\n"
				<< "demog_eta and demog_b (double vector type): " << demog_eta << "\t" << demog_b << "\n"
				<< "Beta values (double vector type): " << TAGInfant << "\t" << TAGPreSAC << "\t" << TAGSAC << "\t" << TAGAdult << "\n"
				<< "Treatment start time, end time and frequency (int vector type): " << treatStart << "\t" << treatEnd << "\t" << treatFreq << "\n"
				<< std::flush;

	/////////////////////////////////////////////////////////////////////////////////////
	//// DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE ///


	// Allocate memory
	dt = (float) 1/nOutputsPerYear;
	nTimeSteps = ((int) ceil(nYears/dt)) + 1;
	results = new wormBurden*[nRepetitions];
	for (int i=0;i<nRepetitions;i++)
	{
		// Allocate repetition
		results[i] = new wormBurden[nTimeSteps];
		memset(results[i],0,sizeof(wormBurden)*nTimeSteps);

		// Initialise nWorms in each repetition
		results[i][0].nWorms = 0; // CHANGE THIS VALUE LATER
		results[i][0].time = 0;
	}


	return true;
}

// Run simulation
void CSimulator::runSimulation()
{
	int runIndex, timeIndex;
	for (runIndex=0;runIndex<nRepetitions;runIndex++)
	{
		wormBurden* currentRun = results[runIndex];
		for (timeIndex=0;timeIndex<nTimeSteps-1;timeIndex++)
		{
			currentRun[timeIndex+1].nWorms = currentRun[timeIndex].nWorms + 1;
			currentRun[timeIndex+1].time = currentRun[timeIndex].time + dt;
		}
	}
}

// Output simulation
void CSimulator::outputSimulation(int n)
{
	for (int i=0;i<nTimeSteps;i++)
	{
		resultsStream << results[n][i].time << "\t"
				<< results[n][i].nWorms << "\n";
	}
	resultsStream << std::flush;
}

//////////////////////////////////////////////////////////////////////////////////
/// Auxiliary function definitions

double* CSimulator::readDoublesVector(char* currentString, int& currentLength)
{
	char* endPointer; // endPointer for each call to strto_ type functions

	// Read in the length of the vector
	int length = strlen(currentString);
	currentLength = length;

	if(length<0)
		return NULL;

	// Create array of doubles
	double* doubleVectorArray = new double[length];
	doubleVectorArray[0] = strtod(currentString,&endPointer);
	for(int i=1;i<length;i++)
	{
		doubleVectorArray[i] = strtod(endPointer, &endPointer);
	}

	return doubleVectorArray;
}

int* CSimulator::readIntsVector(char* currentString, int& currentLength)
{
	char* endPointer; // endPointer for each call to strto_ type functions

	// Read in the length of the vector
	int length = strlen(currentString);
	currentLength = length;

	if(length<0)
		return NULL;

	// Create array of integers
	int* intVectorArray = new int[length];
	intVectorArray[0] = strtol(currentString, &endPointer,10);
	for(int i=1;i<length;i++)
	{
		intVectorArray[i] = strtol(endPointer, &endPointer,10);
	}
	return intVectorArray;
}

// This function takes a random number (0-1) and multiplies it by the max of the passed array.
// It then finds the smallest index that has an array value greater than the product above.
// For a cumulative multinomial array, this will return the index of the event that occurred.
// Also used for drawing a lifespan from the survival curve integral.
int CSimulator::multiNomBasic(double* array, int length, double randNum)
{
	int loopMax = ceil(log(length)/log(2)+2);
	int bottom = -1;
	//double bottomVal;
	int top = length-1;
	double topVal = array[top];
	double target = topVal*randNum;
	int count = 0;

	while(++count<loopMax && (top-bottom>1))
	{
		int mid = (top+bottom)/2;
		double midVal = array[mid];
		if(midVal>=target)
		{
			top = mid;
			topVal = midVal;
		} else
		{
			bottom = mid;
			//bottomVal = midVal;
		}
	}

	if(count>=loopMax)
	{
		logStream << "Max iterations exceeded in multiNomBasic(...),\n" << std::flush;
		return -1;
	}

	return top;
}

int currentRandIndex = 1; // Test line, remove later

// Draw a life span from the survival curve from the population.
double CSimulator::drawLifespan()
{
	// Get a random integer from survivalCurveIntegral using the multinomial generator. This shouldn't be zero!!
	double currentRand = genunf(0,1);

	 // Test lines, remove later
	printf("currentRand number %d: %g\n",currentRandIndex,currentRand);
	currentRandIndex += 1;

	int index = multiNomBasic(survivalCurveIntegral, survivalMaxIndex,currentRand);

	// Interpolate from the returned value.
	double target = currentRand*survivalCurveIntegral[survivalMaxIndex-1];
	double a = (target - survivalCurveIntegral[index-1])/(survivalCurveIntegral[index] - survivalCurveIntegral[index-1]);
	double ans = a*survivalDt + (index - 1)*survivalDt;
	return ans;
}
