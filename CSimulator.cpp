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

	CAGInfant = 0;
	CAGPreSAC = 0;
	CAGSAC = 0;
	CAGAdult = 0;
	InfantBeta = 0;
	PreSACBeta = 0;
	SACBeta = 0;
	AdultBeta = 0;
	R0 = 0;
	lambda = 0;
	gamma = 0;
	z = 0;
	k = 0;
	sigma = 0;
	LDecayRate = 0;
	TAGInfant = 0;
	TAGPreSAC = 0;
	TAGSAC = 0;
	TAGAdult = 0;
	InfantCoverage = 0;
	PreSACCoverage = 0;
	SACCoverage = 0;
	AdultCoverage = 0;
	drugEfficacy = 0;
	treatStart = 0;
	treatEnd = 0;
	treatFreq = 0;
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
	// Number of repetitions
	nRepetitions = atoi(myReader.getParamString("numberRepetitions"));

	// Number of years to run
	nYears = atoi(myReader.getParamString("numberYears"));

	// Number of outputs per year
	nOutputsPerYear = atoi(myReader.getParamString("numberOutputsPerYear"));

	// Number of hosts
	nHosts = atoi(myReader.getParamString("numberHosts"));

	/////////////////////////////////////////////////////////////////////////////
	///  Set up demography

	// Read a parameter to determine whether to use data for the survival curve or a named function

	// Use the exponential-exponential function to define the survival curve
	char* demographyEndPointer;
	demog_eta = strtod(myReader.getParamString("demography"),&demographyEndPointer);
	demog_b = strtod(demographyEndPointer,NULL);
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
		hostPopulation[i] -> birthDate = -myRand()*lifespan;
		hostPopulation[i] -> deathDate = hostPopulation[i] -> birthDate + lifespan;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	//// DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE DEBUG CODE ///

	// Contact age groups
	char * CAGEndPointer;
	CAGInfant = strtol(myReader.getParamString("contactAgeGroups"),&CAGEndPointer,10); // Convert string to integer, base 10
	CAGPreSAC = strtol(CAGEndPointer,&CAGEndPointer,10);
	CAGSAC = strtol(CAGEndPointer,&CAGEndPointer,10);
	CAGAdult = strtol(CAGEndPointer,NULL,10);

	// Beta values
	char * betaEndPointer;
	InfantBeta = strtof(myReader.getParamString("beta"),&betaEndPointer); // Convert string to float
	PreSACBeta = strtof(betaEndPointer,&betaEndPointer);
	SACBeta = strtof(betaEndPointer,&betaEndPointer);
	AdultBeta = strtof(betaEndPointer,NULL);

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
	char * TAGEndPointer;
	TAGInfant = strtol(myReader.getParamString("treatmentAgeGroups"),&TAGEndPointer,10); // Convert string to integer, base 10
	TAGPreSAC = strtol(TAGEndPointer,&TAGEndPointer,10);
	TAGSAC = strtol(TAGEndPointer,&TAGEndPointer,10);
	TAGAdult = strtol(TAGEndPointer,NULL,10);

	// Coverages
	char * coverageEndPointer;
	InfantCoverage = strtof(myReader.getParamString("coverages"),&coverageEndPointer); // Convert string to float
	PreSACCoverage = strtof(coverageEndPointer,&coverageEndPointer);
	SACCoverage = strtof(coverageEndPointer,&coverageEndPointer);
	AdultCoverage = strtof(coverageEndPointer,NULL);

	// Drug efficacy
	drugEfficacy = atof(myReader.getParamString("drugEfficacy"));

	// Chemotherapy timings
	char * chemoTimingsEndPointer;
	treatStart = strtol(myReader.getParamString("treatmentTimes"),&chemoTimingsEndPointer,10); // Convert string to integer, base 10
	treatEnd = strtol(chemoTimingsEndPointer,&chemoTimingsEndPointer,10);
	treatFreq = strtol(chemoTimingsEndPointer,NULL,10);

	// Print some parameters of different types to log file to check if they are being read in correctly
	logStream << "\nRepetitions: " << nRepetitions << "\n"
				<< "Fecundity parameter: " << z << "\n"
				<< "demog_eta and demog_b: " << demog_eta << "\t" << demog_b << "\n"
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
	// Read in the length of the vector
	char* endPointer;
	int length = strtol(currentString,&endPointer,10);
	currentLength = length;

	if(length<0)
		return NULL;

	// Create array of doubles
	double* vectorArray = new double[length];
	for(int i=0;i<length;i++)
	{
		vectorArray[i] = strtod(endPointer, &endPointer);
	}
	return vectorArray;
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

// Uniform random number generator. Should be using the <random> functions.
double CSimulator::myRand()
{
	return (1.0*rand())/RAND_MAX;
}

// Draw a life span from the survival curve from the population.
double CSimulator::drawLifespan()
{
	// Get a random integer from the survivalcurve integral from the multinomial generator. This shouldn't be zero!!
	double currentRand = myRand();
	int index = multiNomBasic(survivalCurveIntegral, survivalMaxIndex,currentRand);

	// Interpolate from the returned value.
	double target = currentRand*survivalCurveIntegral[survivalMaxIndex-1];
	double a = (target - survivalCurveIntegral[index-1])/(survivalCurveIntegral[index] - survivalCurveIntegral[index-1]);
	double ans = a*survivalDt + (index - 1)*survivalDt;
	return ans;
}
