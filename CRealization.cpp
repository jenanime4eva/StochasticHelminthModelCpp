/*
 * File Name: CRealization.cpp
 *
 *  Created on: 14 April 2015
 *  Author: Jie Yang
 *
 */

#include "CRealization.h"
#include "CSimulator.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <vector.h>
#include "randlib.h"

// Class constructor
CRealization::CRealization()
{
	hostPopulation = NULL;
	owner = NULL;
	nHosts = 0;
	tinyIncrement = 0.01;
	sumTotalWorms = 0;
	sumFemaleWorms = 0;
	surveyResultsArray = NULL;
	contactAgeGroupIndex = NULL;
	treatmentAgeGroupIndex = NULL;
}

// Class destructor
CRealization::~CRealization()
{
	int i;
	if(hostPopulation!=NULL)
	{
		for(i=0;i<nHosts;i++)
			delete hostPopulation[i];

		delete[] hostPopulation;
	}

	if (contactAgeGroupIndex != NULL)
		delete[] contactAgeGroupIndex;

	if (treatmentAgeGroupIndex != NULL)
		delete[] treatmentAgeGroupIndex;

	if(surveyResultsArray!=NULL)
	{
		for(i=0;i<owner->surveyResultTimesLength;i++)
			delete[] surveyResultsArray[i];

		delete[] surveyResultsArray;
	}
}

// Set up the realisation
bool CRealization::initialize(CSimulator* currentOwner)
{
	// Set the pointer to the owner
	owner = currentOwner;

	nHosts = owner->nHosts;

	// Set up some arrays
	double* hostTotalWorms = new double[nHosts];
	double* hostFemaleWorms = new double[nHosts];
	double* hostInfectionRate = new double[nHosts];
	double* rates = new double[nHosts+1];

	double sumRates = 0;

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Host lifespans
		double lifespan = owner->drawLifespan();
		hostPopulation[i]->birthDate = -owner->myRandUniform()*lifespan;
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;

		// Add event
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);

		// Total and female worms per individual host
		hostPopulation[i]->totalWorms = owner->myRandPoisson();
		hostPopulation[i]->femaleWorms = owner->myRandBinomial(hostPopulation[i]->totalWorms,0.5); // Half of worm population should be female
		hostTotalWorms[i] = hostPopulation[i]->totalWorms;
		hostFemaleWorms[i] = hostPopulation[i]->femaleWorms;
		sumTotalWorms += hostPopulation[i]->totalWorms;
		sumFemaleWorms += hostPopulation[i]->totalWorms;
	}

	int ratesLength = 1;
	double freelivingNumber = freelivingWorms(hostTotalWorms,hostFemaleWorms,0);
	for (int i=0;i<nHosts;i++)
	{
		hostInfectionRate[i] = freelivingNumber*owner->myRandGamma()*owner->betaValues[hostContactAgeGroupIndex()[i]];
		rates[i] = hostInfectionRate[i];
		sumRates += rates[i];
		ratesLength++; // To obtain the length of the rates array
	}

	// Calculate the event rates (host infection rate and worm total death rate)
	double wormTotalDeathRate = owner->sigma*sumTotalWorms;
	rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array
	sumRates = sumRates + rates[nHosts]; // Sum of the rates

	// Determine if dead worm or new worm

	int event = owner->multiNomBasic2(rates,sumRates,ratesLength,owner->myRandUniform());

	if(event==ratesLength)
	{
		// Dead worm
		int deathIndex = owner->multiNomBasic2(hostTotalWorms,sumTotalWorms,nHosts,owner->myRandUniform());

		// Is this worm female?
		if(owner->myRandUniform()<(hostFemaleWorms[deathIndex]/hostTotalWorms[deathIndex]))
		{
			hostPopulation[deathIndex]->femaleWorms = hostPopulation[deathIndex]->femaleWorms - 1;
		}
		hostPopulation[deathIndex]->totalWorms = hostPopulation[deathIndex]->totalWorms - 1;
	}
	else
	{
		// New worm
		hostPopulation[event]->totalWorms = hostPopulation[event]->totalWorms + 1;
		printf("\nhostPopulation[event]->totalWorms: %f",hostPopulation[event]->totalWorms);
		// Female worm
		if(owner->myRandUniform()<0.5)
		{
			hostPopulation[event]->femaleWorms = hostPopulation[event]->femaleWorms + 1;
		}
	}

	// Delete any arrays not in use anymore
	if (hostTotalWorms != NULL)
		delete[] hostTotalWorms;

	if (hostFemaleWorms != NULL)
		delete[] hostFemaleWorms;

	if (hostInfectionRate != NULL)
		delete[] hostInfectionRate;

	if (rates != NULL)
		delete[] rates;

	// USAGE: addEvent(int newType, double newTime, void* currentSubject = NULL)


	// Add treatment events
	double treatmentTimeElapsed = 0;
	while((owner->treatEnd - owner->treatStart) != treatmentTimeElapsed)
	{
		localEvents.addEvent(CHEMOTHERAPY,owner->treatStart+treatmentTimeElapsed,NULL);
		treatmentTimeElapsed = treatmentTimeElapsed + owner->treatInterval;
	}
	localEvents.addEvent(CHEMOTHERAPY,owner->treatEnd,NULL); // Last treatment time

	// Add run termination point
	localEvents.addEvent(TERMINATE,owner->nYears,NULL);

	// Set up results collection for this realisation
	if(owner->surveyResultTimes!=NULL)
	{
		// There are survey results to collect
		// For each time point set up an array of surveyResultData structures the length of the population size
		surveyResultsArray = new surveyResultData*[owner->surveyResultTimesLength];
		for(int i=0;i<owner->surveyResultTimesLength;i++)
		{
			surveyResultsArray[i] = new surveyResultData[nHosts];
			// Add a pointer to the array that's going to hold the data
			localEvents.addEvent(SURVEY_EVENT,owner->surveyResultTimes[i],surveyResultsArray[i]);
		}
	}

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,owner->nYears - 5.0,NULL);

	return true;
}

// Run the realisation
bool CRealization::run()
{
	Event currentEvent;

	do
	{
		// Get next predetermined event
		localEvents.popEvent(currentEvent);

		// Do stochastic events up to the time of the next event

		// Do the event
		switch (currentEvent.type)
		{
			case HOST_DEATH:
				hostDeathResponse(currentEvent);
				break;
			case CHEMOTHERAPY:
				treatmentResponse(currentEvent);
				printf("\nCHEMOTHERAPY at time: %f",currentEvent.time); // Print test flag to console
				break;
			case DEBUG_EVENT:
				debugEventResponse(currentEvent);
				break;
			case SURVEY_EVENT:
				surveyResultResponse(currentEvent);
				break;
			case TERMINATE:
				break;
			default:
				owner->logStream << "Event number " << currentEvent.type << " not known.\n"
				<< "Look at #define values in CPreDetEventQueue.h for definition.\n" << std::flush;
				break;
		}
	} while(currentEvent.type!=TERMINATE); // Do events until the last survey year is reached

	return true;
}

// Respond to host death
bool CRealization::hostDeathResponse(Event& currentEvent)
{
	// Rejuvenate host

	CHost* currentHost = (CHost*) currentEvent.subject;

	// Set birth date to now
	currentHost->birthDate = currentEvent.time;

	// Calculate death dates
	double lifespan = owner->drawLifespan();
	currentHost->deathDate = currentEvent.time + lifespan;

	// Kill off all their worms
	currentHost->totalWorms = 0;
	currentHost->femaleWorms = 0;

	// Assign death event
	localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

	return true;
}

// Respond to treatment
bool CRealization::treatmentResponse(Event& currentEvent)
{
	CHost* currentHost = (CHost*) currentEvent.subject;

	int currentContactIndex = 0;
	int individualTreatedIndex = 0;

	// Contact age group index for current host
	while (owner->contactAgeBreaks[currentContactIndex] + tinyIncrement < -currentHost->birthDate)
	{
			currentContactIndex++;
	}
	currentContactIndex = currentContactIndex-1;

	double individualCoverage = owner->coverage[currentContactIndex];

	// Index for whether current host is treated, allocate value of 1 if statement is fulfilled else -1
	if(owner->myRandUniform()<individualCoverage)
	{
		individualTreatedIndex = 1;
	}
	else
	{
		individualTreatedIndex = -1;
	}

	// If individualTreatedIndex is NOT equal to -1
	if (individualTreatedIndex != -1)
	{
		// How many worms to die?
		double totalW = currentHost->totalWorms;
		double femaleW = currentHost->femaleWorms;
		double maleW = totalW - femaleW;

		double maleToDie = owner->myRandBinomial(maleW,owner->drugEff);
		double femaleToDie = owner->myRandBinomial(femaleW,owner->drugEff);

		currentHost->totalWorms = totalW - maleToDie - femaleToDie;
		currentHost->femaleWorms = femaleW - femaleToDie;
	}

	return true;
}

// What to output from the simulation
bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// Collect data from each host individual
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	/*
	// For looking at the host ages across time
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].age = currentEvent.time - hostPopulation[i]->birthDate;
	}
	*/

	// For looking at mean total worm burdens for individual runs across time
	double sumTotalWormsPerRun = 0;
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].nTotalWorms = hostPopulation[i]->totalWorms;
		sumTotalWormsPerRun = outputArray[i].nTotalWorms + sumTotalWormsPerRun; // Get sum of total worms per run
	}
	outputArray->meanTotalWormsPerRun = sumTotalWormsPerRun/nHosts;

	// For looking at mean female worm burdens for individual runs across time
	double sumFemaleWormsPerRun = 0;
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].nFemaleWorms += hostPopulation[i]->femaleWorms;
		sumFemaleWormsPerRun = outputArray[i].nFemaleWorms + sumFemaleWormsPerRun; // Get sum of female worms per run
	}
	outputArray->meanFemaleWormsPerRun = sumFemaleWormsPerRun/nHosts;

	return true;
}

// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{

}


//////////////////////////////////////////////////////////////////////////////////
/// OTHER FUNCTIONS

// ts = time do the event at
double CRealization::freelivingWorms(double* totalWormArray, double* femaleWormArray, double ts)
{
	double freeliving = 4; // Don't know what this number should be
	double sumEggsOutputPerHostRhoValues = 0;
	double eggsProductionRate = 0;
	double* productiveFemaleWorms = new double[nHosts];
	double* eggsOutputPerHost = new double[nHosts];

	for(int i=0;i<nHosts;i++)
	{
		productiveFemaleWorms[i] = femaleWormArray[i];
		// Female worms produce fertilised eggs only if there is a male worm around
		if (totalWormArray[i]==femaleWormArray[i])
		{
			productiveFemaleWorms[i] = 0;
		}

		eggsOutputPerHost[i] = owner->lambda*productiveFemaleWorms[i]*exp(-productiveFemaleWorms[i]*owner->gamma);

		sumEggsOutputPerHostRhoValues += eggsOutputPerHost[i]*owner->rhoValues[hostContactAgeGroupIndex()[i]];
	}

	eggsProductionRate = owner->psi*sumEggsOutputPerHostRhoValues/(owner->nHosts);

	double expo = exp(-owner->ReservoirDecayRate*ts);

	freeliving = freeliving*expo + eggsProductionRate*(1-expo)/(owner->ReservoirDecayRate);

	// Delete stuff we don't need anymore
	if (productiveFemaleWorms != NULL)
		delete[] productiveFemaleWorms;

	if (eggsOutputPerHost != NULL)
		delete[] eggsOutputPerHost;

	return freeliving;
}

int* CRealization::hostContactAgeGroupIndex()
{
	int currentContactIndex = 0;
	contactAgeGroupIndex = new int[nHosts];

	// Contact age group index for each host
	for(int i=0;i<nHosts;i++)
	{
		while (owner->contactAgeBreaks[currentContactIndex] + tinyIncrement < -hostPopulation[i]->birthDate)
		{
			currentContactIndex++;
		}
		contactAgeGroupIndex[i] = currentContactIndex-1;
	}

	return contactAgeGroupIndex;
}

int* CRealization::hostTreatmentAgeGroupIndex()
{
	int currentTreatmentIndex = 0;
	treatmentAgeGroupIndex = new int[nHosts];

	// Treatment age group index for each host
	for(int i=0;i<nHosts;i++)
	{
		while (owner->treatmentBreaks[currentTreatmentIndex] + tinyIncrement < -hostPopulation[i]->birthDate)
		{
			currentTreatmentIndex++;
		}
		treatmentAgeGroupIndex[i] = currentTreatmentIndex-1;
	}

	return treatmentAgeGroupIndex;
}
