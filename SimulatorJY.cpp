 #include "SimulatorJY.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "randlib.h"

helsimJY::helsimJY(void)
{
	results = NULL; 
}

helsimJY::~helsimJY(void)
{
	// Delete results allocations, last thing.
	if(results!=NULL)
	{
		// Need to delete memory...
		for(int i=0;i<nRepetitions;i++)
		{
			delete[] results[i]; 
		}

		delete[] results; 
	}
}


// Initialise the input/output aspects of the simulator
bool helsimJY::initialiseIO(char* logFileName, char* paramFileName, char* resultsFileName)
{
	// Setting up all the files we'll need
	logStream.open(logFileName); 
	if(!logStream.is_open())
	{
		std::cout << "Couldn't open the log file: " << logFileName << "\nexiting\n"; 
		return false; 
	}

	paramsIn.open(paramFileName); 
	if(!paramsIn.is_open())
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

	// Everything's ok, apparently
	return true; 
}

// General initialisation
bool helsimJY::initialiseSimulation()
{
	// Read in parameters line by line

	// SIMULATION STUFF
	if(!getNextLine(&paramsIn))
		return false; 
	nRepetitions = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false; 
	nMaxYear = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	nOutFreq = atoi(buffer);

	// HOST STUFF
	if(!getNextLine(&paramsIn))
		return false;
	hostPopulation = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false; 
	demogEta = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	demogB = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	nAgeGroups = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	contactBreakInfants = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	contactBreakPreSAC = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	contactBreakSAC = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	contactBreakAdults = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	InfantBeta = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	PreSACBeta = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	SACBeta = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	AdultBeta = atof(buffer);

	// WORM STUFF
	if(!getNextLine(&paramsIn))
		return false;
	R0 = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	lambda = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	gamma = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	z = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	k = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	sigma = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	LDecayRate = atof(buffer);

	// TREATMENT STUFF
	if(!getNextLine(&paramsIn))
		return false;
	treatBreakInfants = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	treatBreakPreSAC = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	treatBreakSAC = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	treatBreakAdults = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	InfantCoverage = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	PreSACCoverage = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	SACCoverage = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	AdultCoverage = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	drugEfficacy = atof(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	treatStartYear = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	treatEndYear = atoi(buffer);

	if(!getNextLine(&paramsIn))
		return false;
	treatFreq = atoi(buffer);

	// Set up the system
	// Allocate memory...
	nTimeSteps = ((int) ceil(nMaxYear/tstep)) + 1; // TO DO: DEFINE tstep
	results = new meanWormBurden*[nRepetitions];
	for(int i=0;i<nRepetitions;i++)
	{
		// Allocate results array for each repetition
		results[i] = new meanWormBurden[nTimesteps];
		memset(results[i], 0, sizeof(meanWormBurden)*nTimesteps);

		// Initialise for each repetition
		results[i][0].M = initM; // Initial mean worm burden
		results[i][0].time = 0; 
	}

	return true; 
}

// Run the simulation
void helsimJY::runSimulation()
{
	int nInfections, nRecoveries; 
	int rIndex, timeIndex; 
	for(rIndex=0;rIndex<nRepetitions;rIndex++)
	{
		meanWormBurden* currentRun = results[rIndex];
		for(timeIndex=0;timeIndex<nTimesteps-1;timeIndex++)
		{
			nInfections = nRecoveries = 0;
			if(currentRun[timeIndex].inf!=0)
			{
				nInfections = ignbin(currentRun[timeIndex].susc,beta*currentRun[timeIndex].inf*dt/totalPop); 
				nRecoveries = ignbin(currentRun[timeIndex].inf, sigma*dt); 
			}

			currentRun[timeIndex+1].susc = currentRun[timeIndex].susc - nInfections; 
			currentRun[timeIndex+1].inf = currentRun[timeIndex].inf + nInfections - nRecoveries; 
			currentRun[timeIndex+1].rec = currentRun[timeIndex].rec + nRecoveries; 
			currentRun[timeIndex+1].time = currentRun[timeIndex].time + dt; 
		}
	}
}

// Output repetition number n
void CSimulator::outputRealisation(int n)
{
	for (int i=0;i<nTimesteps;i++)
	{
		resultsStream << results[n][i].time << "\t"
			<< results[n][i].susc << "\t"
			<< results[n][i].inf << "\t"
			<< results[n][i].rec << "\n";
	}
}

void CSimulator::outputPropExtinct()
{
	for(int i=0;i<nTimesteps;i++)
	{
		double propExtinct = 0; 
		for(int j=0;j<nRealisations;j++)
		{
			if(results[j][i].inf==0)
				propExtinct +=1.0; 
		}
		propExtinct /= nRealisations; 

		resultsStream << results[0][i].time << "\t" << propExtinct << "\n"; 
	}
}

// Get the next usable line from this stream
bool CSimulator::getNextLine(std::ifstream* currentStream)
{
	// I'm passing a pointer to the stream in case I ever need to read from some other file that isn't param file. 
	if(currentStream->eof())
		return false; // Already at the end, nothing more here...

	currentStream->getline(buffer, BUFFER_SIZE); 

	// If NOT at the end of file AND line is rubbish, get another line...
	while(!currentStream->eof()&&isLineRubbish(buffer))
	{
		currentStream->getline(buffer,BUFFER_SIZE);
	}

	if(isLineRubbish(buffer))
		return false; // Must have exited while loop due to end of file, hence run out of lines

	return true; 
}

// Is this line rubbish?
bool CSimulator::isLineRubbish(char* currentLine)
{
	// This line is rubbish if length = 0 OR first character = %
	return (strlen(currentLine)==0)||(currentLine[0]=='%');
}
