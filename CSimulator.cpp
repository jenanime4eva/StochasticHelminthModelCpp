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
//#include <string>


// Class Constructor
CSimulator::CSimulator()
{
	// Default values
	nRepetitions = 0;
	nYears = 0;
	nOutputsPerYear = 0;
	nTimeSteps = 0;
	dt = 0;
	nHosts = 0;
	startYear = 0;

	survivalCurve = survivalCurveCumul = NULL;
	hostMu = probDeath = probDeathIntegral = NULL;
	demogDt = 0;
	demog_b = demog_eta = 0;
	hostMuData = NULL;
	hostMuDataLength = 0;
	muDataUpperBounds = NULL;
	muUpperBoundsLength = 0;
	upperAgeBound = 0;
	maxDtIntevals = 0;

	// Results-related variables.
	surveyResultTimes = NULL;
	surveyResultTimesLength = 0;



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
	if(survivalCurve!=NULL)
		delete[] survivalCurve;

	if(survivalCurveCumul!=NULL)
		delete[] survivalCurveCumul;

	if(surveyResultTimes!=NULL)
		delete[] surveyResultTimes;

	if(hostMuData!=NULL)
		delete[] hostMuData;

	if(muDataUpperBounds!=NULL)
		delete[] muDataUpperBounds;

	if(hostMu!=NULL)
		delete[] hostMu;

	if(probDeath!=NULL)
		delete[] probDeath;

	if(probDeathIntegral!=NULL)
		delete[] probDeathIntegral;
}

// Initialise the input/output aspects of the simulator
bool CSimulator::initialiseIO(char* run, char* path, char* paramFilePath)
{
	// as strings classes...
	runName = run;
	thePath = path;

	// Setting up all the files needed.
	std::string logFilePath = thePath + runName + ".log.txt";
	logStream.open(logFilePath.c_str());
	if(!logStream.is_open())
	{
		std::cout << "Couldn't open the log file: " << logFilePath << "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Log file opened.\n" << std::flush;

	myReader.setNewFileName(paramFilePath);
	/*
	if (!myReader.setNewFileName(paramFilePath))
 	{
		logStream << "Couldn't open the param file: " << paramFilePath << "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Param reader configured.\n" << std::flush;
	*/

	// stub name for all results files.
	resultsStub = thePath + runName;


	// Everything is ok
	return true;
}

// General initialisation
bool CSimulator::initialiseSimulation()
{
	char* endPointer; // general variable for the end pointer used in strto_ functions.
	char* temp; // general purpose pointer to string.
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

	// use the expo-expo function to define the survival curve.   DEBUGDEBUG some of these may not be needed anymore.
	demog_eta = strtod(myReader.getParamString("demog_eta"),&endPointer);
	demog_b = strtod(myReader.getParamString("demog_b"),&endPointer);
	demogDt = strtod(myReader.getParamString("demogDt"),&endPointer);


	////////////////////////////////////////////////////////////////////////////////

	// read in death rates.
	temp = myReader.getParamString("hostMu");
	if(temp!=NULL)
	{
		hostMuData = readDoublesVector(temp,hostMuDataLength);
	}

	// read in death rate upper bounds.
	temp = myReader.getParamString("muUpperB");
	if(temp!=NULL)
	{
		muDataUpperBounds = readDoublesVector(temp, muUpperBoundsLength);
	}


	// construct the mu, prob of death and survival vectors.
	//double upperAgeBound = muUpperBounds[muUpperBoundsLength-1];
	maxDtIntevals = (int) floor(muDataUpperBounds[muUpperBoundsLength-1]/demogDt);
	upperAgeBound = maxDtIntevals*demogDt;
	int currentMuIndex = 0;
	double currentSurvival = 1;
	double currentSurvivalCumul = 0;
	double currentProbDeathCumul = 0;
	double currentMuDtCumul = 0;
	double tiny = 0.01;
	survivalCurve = new double[maxDtIntevals];
	survivalCurveCumul = new double[maxDtIntevals];
	hostMu = new double[maxDtIntevals];
	probDeath = new double[maxDtIntevals];
	probDeathIntegral = new double[maxDtIntevals];

	for(int i=0;i<maxDtIntevals;i++)
	{
		double currentIntEnd = (i+1)*demogDt;
		if(muDataUpperBounds[currentMuIndex]+tiny<currentIntEnd) // is current dt interval within data upper bound?
			currentMuIndex++;
		hostMu[i] = hostMuData[currentMuIndex];
		probDeath[i] = currentSurvival*hostMu[i]*demogDt;
		currentMuDtCumul += hostMu[i]*demogDt;

		probDeathIntegral[i] = currentProbDeathCumul + probDeath[i];
		currentProbDeathCumul = probDeathIntegral[i];

		survivalCurve[i] = exp(-currentMuDtCumul);
		currentSurvival = survivalCurve[i];

		survivalCurveCumul[i] = survivalCurve[i] + currentSurvivalCumul;
		currentSurvivalCumul = survivalCurveCumul[i];
	}


	//////////// Set up results collection.
	temp = myReader.getParamString("surveyTimes");
	if(temp!=NULL)
	{
		// there are survey results to collect.
		surveyResultTimes = readDoublesVector(temp,surveyResultTimesLength);
	}

	//////////// Test area ////////////////////////////////////////
	// set up a realization.
	myRealization.initialize(this);
	myRealization.run();


	return true;
}

// Run simulation
void CSimulator::runSimulation()
{
}

// Output simulation
void CSimulator::outputSimulation()
{
	// output survey-type results if there's any...
	if(surveyResultTimes!=NULL)
	{
		// there's some survey results...
		// create output stream...
		std::string surveyResultsOut = thePath + runName + ".surveyResults.txt";
		std::ofstream surveyStream(surveyResultsOut.c_str());

		// DEBUG: currently only ONE realization, so...
		for(int j=0;j<surveyResultTimesLength;j++)
			surveyStream << surveyResultTimes[j] << "\t";
		surveyStream << "\n";
		// loop through the hosts...
		for(int i=0;i<nHosts;i++)
		{
			for(int j=0;j<surveyResultTimesLength;j++)
			{
				surveyStream << myRealization.surveyResultsArray[j][i].age << "\t";
			}
			surveyStream << "\n";
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
/// Auxiliary function definitions.

// creates a vector of doubles from string from param file.
// Allocates memory, so need to delete.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TO DO: Write this code so it doesn't need to know vector length, to match the format for the R code.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double* CSimulator::readDoublesVector(char* currentString, int& currentLength)
{
	// read in the length of the vector.
	char* endPointer;

	//// Count the length of the vector. It's important for the string to not have trailing zeros.
	// These are currently removed in the param reading function.
	int count = 1;
	double junk = strtod(currentString, &endPointer);
	while(strlen(endPointer)>0)
	{
		junk = strtod(endPointer, &endPointer);
		count++;
	}

	currentLength = count;

	// create array of doubles.
	double* vectorArray = new double[currentLength];
	vectorArray[0] = strtod(currentString, &endPointer);
	for(int i=1;i<currentLength;i++)
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
	// get a random integer from the probDeathIntegral integral from the multinomial generator. This shouldn't be zero!!
	double currentRand = myRand();
	int index = multiNomBasic(probDeathIntegral, maxDtIntevals,currentRand);

	double ans = (index+0.5)*demogDt;  // choose a point in the middle of the interval in which the person dies.
	return ans;
}
