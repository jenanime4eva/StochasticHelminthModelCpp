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

// Class Constructor
CSimulator::CSimulator()
{
	// Default values
	hostPopulation = NULL;
	results = NULL;
	nRepetitions = 0;
	nYears = 0;
	nOutputsPerYear = 0;
	nHosts = 0;
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
	nTimeSteps = 0;
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
	logStream << "Log file opened.\n";

	if (!myReader.setNewFileName(paramFileName))
 	{
		logStream << "Couldn't open the param file: " << paramFileName << "\nexiting\n";
		return false;
	}

	resultsStream.open(resultsFileName);
	if(!resultsStream.is_open())
	{
		logStream << "Couldn't open the results file: " << resultsFileName << "\nexiting\n";
		return false;
	}

	// Get IO related parameters
	nRepetitions = atoi(myReader.getParamString("param1"));
	nYears = atoi(myReader.getParamString("param2"));
	nOutputsPerYear = atoi(myReader.getParamString("param3"));

	//Test: Parameter printout to console
	printf ("Repetitions: %d\n",nRepetitions);
	printf ("Years to run: %d\n",nYears);
	printf ("Outputs per year: %d\n",nOutputsPerYear);


	// Everything is ok
	return true;
}

// General initialisation
bool CSimulator::initialiseSimulation()
{
	// Get model related parameters

	// Number of hosts
	nHosts = atoi(myReader.getParamString("param4"));

	// Demography info
	char * demog;
	demog_eta = strtof(myReader.getParamString("param5"),&demog); // Convert string to float
	demog_b = strtof(demog,NULL);

	// Number of age groups
	nAG = atoi(myReader.getParamString("param6"));

	// Contact age groups
	char * CAGEndPointer;
	CAGInfant = strtol(myReader.getParamString("param7"),&CAGEndPointer,10); // Convert string to integer, base 10
	CAGPreSAC = strtol(CAGEndPointer,&CAGEndPointer,10);
	CAGSAC = strtol(CAGEndPointer,&CAGEndPointer,10);
	CAGAdult = strtol(CAGEndPointer,NULL,10);

	// Beta values
	char * betaEndPointer;
	InfantBeta = strtof(myReader.getParamString("param8"),&betaEndPointer); // Convert string to float
	PreSACBeta = strtof(betaEndPointer,&betaEndPointer);
	SACBeta = strtof(betaEndPointer,&betaEndPointer);
	AdultBeta = strtof(betaEndPointer,NULL);

	// R0
	R0 = atoi(myReader.getParamString("param9"));

	// Egg per gram
	lambda = atoi(myReader.getParamString("param10"));

	// Exponential density dependence of parasite adult stage
	gamma = atof(myReader.getParamString("param11"));

	// Fecundity parameter z = exp(-gamma)
	z = atof(myReader.getParamString("param12"));

	// Shape parameter of assumed negative binomial distribution of worms amongst host
	k = atof(myReader.getParamString("param13"));

	// Worm death rate i.e. 1/worm_life_span
	sigma = atof(myReader.getParamString("param14"));

	// Reservoir decay rate (decay rate of eggs in the environment)
	LDecayRate = atoi(myReader.getParamString("param15"));

	// Treatment age groups
	char * TAGEndPointer;
	TAGInfant = strtol(myReader.getParamString("param16"),&TAGEndPointer,10); // Convert string to integer, base 10
	TAGPreSAC = strtol(TAGEndPointer,&TAGEndPointer,10);
	TAGSAC = strtol(TAGEndPointer,&TAGEndPointer,10);
	TAGAdult = strtol(TAGEndPointer,NULL,10);

	// Coverages
	char * coverageEndPointer;
	InfantCoverage = strtof(myReader.getParamString("param17"),&coverageEndPointer); // Convert string to float
	PreSACCoverage = strtof(coverageEndPointer,&coverageEndPointer);
	SACCoverage = strtof(coverageEndPointer,&coverageEndPointer);
	AdultCoverage = strtof(coverageEndPointer,NULL);

	// Drug efficacy
	drugEfficacy = atof(myReader.getParamString("param18"));

	// Chemotherapy timings
	char * chemoTimingsEndPointer;
	treatStart = strtol(myReader.getParamString("param19"),&chemoTimingsEndPointer,10); // Convert string to integer, base 10
	treatEnd = strtol(chemoTimingsEndPointer,&chemoTimingsEndPointer,10);
	treatFreq = strtol(chemoTimingsEndPointer,NULL,10);

	//Test: Parameter printout to console
	printf ("Host population number: %d\n",nHosts);
	printf ("Minimum age of treatment groups for infants, pre-SACs, SACs and adults : %d %d %d %d\n",TAGInfant,TAGPreSAC,TAGSAC,TAGAdult);
	printf ("Coverages for infants, pre-SAC, SAC and adults: %g %g %g %g\n",InfantCoverage,PreSACCoverage,SACCoverage,AdultCoverage);
	printf ("Drug Efficacy: %g\n",drugEfficacy);
	printf ("Treatment start year, end year and frequency: %d %d %d\n",treatStart,treatEnd,treatFreq);

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;
	}

	// Allocate memory
	nTimeSteps = (int) nYears/nOutputsPerYear;
	results = new wormBurden*[nRepetitions];
	for (int i=0;i<nRepetitions;i++)
	{
		// Allocate repetition
		results[i] = new wormBurden[nTimeSteps];
		memset (results[i],0,sizeof(wormBurden)*nTimeSteps); // LOOK INTO THIS

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
		for (timeIndex=0;timeIndex<nTimeSteps;timeIndex++)
		{
			//currentRun[timeIndex+1].nWorms = 0; // THIS LINE CAUSES CRASH, LOOK INTO IT
			currentRun[timeIndex+1].time = currentRun[timeIndex].time + nOutputsPerYear;
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
}



