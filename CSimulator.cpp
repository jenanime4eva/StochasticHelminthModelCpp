/*
 * CSimulator.cpp
 *
 *  Created on: 26 Nov 2014
 *      Author: Jie Yang
 */

#include ".\CSimulator.h"
#include "CHost.h"
#include "CParamReader.h"
#include "CRealization.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>


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

	survivalCurve = survivalCurveIntegral = NULL;
	survivalDt = 0;
	survivalMaxIndex = 0;
	demog_b = demog_eta = 0;
	/*
	demog_eta = 0;
	demog_b = 0;
	nAG = 0;
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
	*/
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
	char* endPointer; // general variable for the end pointer used in strto_ functions.
	// Get model related parameters
	// Number of repetitions
	nRepetitions = atoi(myReader.getParamString("param1"));

	// when to start (may do this differently later).
	startYear = atof(myReader.getParamString("startYear"));

	// Number of years to run
	nYears = atof(myReader.getParamString("nYears"));

	// Number of hosts
	nHosts = atoi(myReader.getParamString("nHosts"));

	/////////////////////////////////////////////////////////////////////////////
	///  Set up demography.

	// use the expo-expo function to define the survival curve.
	demog_eta = strtod(myReader.getParamString("demog_eta"),&endPointer);
	demog_b = strtod(myReader.getParamString("demog_b"),&endPointer);
	survivalDt = strtod(myReader.getParamString("survivalDt"),&endPointer);
	double survivalMaxAge = strtod(myReader.getParamString("survivalMaxAge"),&endPointer);

	survivalMaxIndex = (int) ceil(survivalMaxAge/survivalDt);
	survivalCurve = new double[survivalMaxIndex];
	survivalCurveIntegral = new double[survivalMaxIndex];
	double subTotal = 0;
	for(int i=0;i<survivalMaxIndex;i++)
	{
		//survivalCurve[i] = exp(-demog_eta*(exp(demog_b*survivalDt*i)-1));
		survivalCurve[i] = exp(-demog_eta*survivalDt*i);
		subTotal += survivalCurve[i];
		survivalCurveIntegral[i] = (subTotal - (survivalCurve[0]+survivalCurve[i])/2)*survivalDt;
	}


	//////////// Test area ////////////////////////////////////////
	// set up a realization.
	CRealization myRealization(this);
	myRealization.initialize();
	myRealization.run();

	/*
	////////////////////////////////// THIS NEEDS TO GO IN THE REALIZATION ////////////////////////////////////
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
	*/

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
			currentRun[timeIndex+1].nWorms = currentRun[timeIndex].nWorms + 1; // THIS LINE CAUSES CRASH, LOOK INTO IT
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
/// Auxiliary function definitions.

double* CSimulator::readDoublesVector(char* currentString, int& currentLength)
{
	// read in the length of the vector.
	char* endPointer;
	int length = strtol(currentString,&endPointer,10);
	currentLength = length;

	if(length<0)
		return NULL;

	// create array of doubles.
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

// uniform random number generator. Should be using the <random> functions.
double CSimulator::myRand()
{
	return (1.0*rand())/RAND_MAX;
}

// draw a life span from the survival curve from he population.
double CSimulator::drawLifespan()
{
	// get a random integer from the survivalcurve integral from the multinomial generator. This shouldn't be zero!!
	double currentRand = myRand();
	int index = multiNomBasic(survivalCurveIntegral, survivalMaxIndex,currentRand);

	// interpolate from the returned value.
	double target = currentRand*survivalCurveIntegral[survivalMaxIndex-1];
	double a = (target - survivalCurveIntegral[index-1])/(survivalCurveIntegral[index] - survivalCurveIntegral[index-1]);
	double ans = a*survivalDt + (index - 1)*survivalDt;
	return ans;
}
