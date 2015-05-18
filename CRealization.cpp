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
	hostContactAgeGroupIndex = NULL;
	hostTreatmentAgeGroupIndex = NULL;
	productiveFemaleWorms = NULL;
	freelivingWorms = NULL;
	hostInfectionRate = NULL;
	rates = NULL;
	individualCoverage = NULL;
	surveyResultsArray = NULL;
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

	if (hostContactAgeGroupIndex != NULL)
		delete[] hostContactAgeGroupIndex;

	if (hostTreatmentAgeGroupIndex != NULL)
		delete[] hostTreatmentAgeGroupIndex;

	if (productiveFemaleWorms != NULL)
		delete[] productiveFemaleWorms;

	if (freelivingWorms != NULL)
		delete[] freelivingWorms;

	if (hostInfectionRate != NULL)
		delete[] hostInfectionRate;

	if (rates != NULL)
		delete[] rates;

	if (individualCoverage != NULL)
		delete[] individualCoverage;

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
	hostContactAgeGroupIndex = new double[nHosts];
	hostTreatmentAgeGroupIndex = new double[nHosts];
	productiveFemaleWorms = new double[nHosts];
	freelivingWorms = new double[nHosts];
	hostInfectionRate = new double[nHosts];
	rates = new double[nHosts];
	individualCoverage = new double[nHosts];

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Host lifespans
		double lifespan = owner->drawLifespan();
		hostPopulation[i]->birthDate = -owner->myRandUniform()*lifespan;
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;

		int currentContactIndex = 0;
		int currentTreatmentIndex = 0;

		// Contact age group index for each host
		while (owner->contactAgeBreaks[currentContactIndex] + tinyIncrement < -hostPopulation[i]->birthDate)
		{
			currentContactIndex++;
		}
		hostContactAgeGroupIndex[i] = currentContactIndex-1;

		// Treatment age group index for each host
		while (owner->treatmentBreaks[currentTreatmentIndex] + tinyIncrement < -hostPopulation[i]->birthDate)
		{
			currentTreatmentIndex++;
		}
		hostTreatmentAgeGroupIndex[i] = currentTreatmentIndex-1;

		// Add event
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);

		// Total and female worms per individual host
		hostPopulation[i]->totalWorms = owner->myRandPoisson();
		hostPopulation[i]->femaleWorms = owner->myRandBinomial();
		sumTotalWorms += hostPopulation[i]->totalWorms;
		sumFemaleWorms += hostPopulation[i]->totalWorms;

		//double freeliving = 4; // Don't know what this number should be

		// Female worms produce fertilised eggs only if there is a male worm around
		//productiveFemaleWorms[i] = hostPopulation[i]->femaleWorms;
		//if (hostPopulation[i]->totalWorms==hostPopulation[i]->femaleWorms)
		//{
			//productiveFemaleWorms[i] = 0;
		//}

		//double eggsOutputPerHost = owner->lambda*productiveFemaleWorms*exp(-productiveFemaleWorms*owner->gamma);

		//use currentEvent.time  for current time, WHAT SHALL TO DO ABOUT TS?


		// Calculate the event rates (host infection rate and worm total death rate)
		//hostInfectionRate[i] = freelivingWorms[i]*owner->myRandGamma()*owner->betaValues[hostContactAgeGroupIndex[i]];
		//hostInfectionRate[i] = freeliving*owner->myRandGamma()*owner->betaValues[hostContactAgeGroupIndex[i]];
		//rates[i] = hostInfectionRate[i];
	}

	//double wormTotalDeathRate = owner->sigma*sumTotalWorms;
	//rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array

	// Determine if dead worm or new worm

	//int ratesLength = sizeof(rates)/sizeof(rates[0]);
	//double randNum = owner->myRandUniform();
	/*
	int event = owner->which(rates,randNum);

	if(event==ratesLength)
	{
		// Dead worm
		//int deathIndex = owner->which(currentHost->totalWorms,randNum);

		// Is this worm female?

	}
	else
	{
		// New worm
	}
	*/




	// Add treatment start point
	localEvents.addEvent(TREATMENT_START,owner->treatStart,NULL);

	/*
	// Coverage levels for each host
	for(int i=0;i<nHosts;i++)
	{
		individualCoverage[i] = owner->coverage[hostContactAgeGroupIndex[i]];
	}
	*/

	// Add treatment end point
	localEvents.addEvent(TREATMENT_END,owner->treatEnd,NULL);


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
			case TREATMENT_START:
				break;
			case TREATMENT_END:
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
				owner->logStream << "Event number " << currentEvent.type << " not known.\n" << std::flush;
				break;
		}
	}while(currentEvent.type!=TERMINATE);

	return true;
}

// Respond to host death
bool CRealization::hostDeathResponse(Event& currentEvent)
{
	// Rejuvenate host
	CHost* currentHost = (CHost*) currentEvent.subject;
	// Set birth date to now
	currentHost->birthDate = currentEvent.time;
	double lifespan = owner->drawLifespan();
	// Calculate death dates
	currentHost->deathDate = currentEvent.time + lifespan;

	// Kill off all their worms
	currentHost->totalWorms = 0;
	currentHost->femaleWorms = 0;

	// Assign death event
	localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

	// Add totalWorms and femaleWorms events here?

	return true;
}

// What to output from the simulation
bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// Collect data from each host individual
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	// For looking at the host ages across time
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].age = currentEvent.time - hostPopulation[i]->birthDate;
	}

	// For looking at total worm burdens for individual hosts across time
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].nTotalWorms = hostPopulation[i]->totalWorms;
	}

	// For looking at female worm burdens for individual hosts across time
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].nFemaleWorms = hostPopulation[i]->femaleWorms;
	}

	return true;
}

// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{
}
