/*
 * File Name: CRealization.cpp
 *
 *  Created on: 14 April 2015
 *  Author: Jie Yang
 *
 */

#include "CRealization.h"
#include "CSimulator.h"

// Class constructor
CRealization::CRealization()
{
	hostPopulation = NULL;
	owner = NULL;
	nHosts = 0;
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

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Host lifespans
		double lifespan = owner->drawLifespan();
		hostPopulation[i]->birthDate = -owner->myRandUniform()*lifespan;
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);

		// Total and female worms per individual host
		hostPopulation[i]->totalWorms = owner->myRandPoisson();
		hostPopulation[i]->femaleWorms = owner->myRandBinomial();
	}

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
			case RATE_EVENT:
				calculateEventRatesResponse(currentEvent);
				break;
			case HOST_DEATH:
				hostDeathResponse(currentEvent);
				break;
			case WORM_BIRTH_DEATH:
				wormBirthDeathResponse(currentEvent);
				break;
			case WORM_FREELIVING:
				wormFreelivingResponse(currentEvent);
				break;
			case TREATMENT_EVENT:
				chemotherapyResponse(currentEvent);
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

// Calculate the event rates (worm total death rate and host infection rate)
bool CRealization::calculateEventRatesResponse(Event& currentEvent)
{

	CHost* currentHost = (CHost*) currentEvent.subject;

	// INSERT DEATHRATE AND HOSTINFECTIONRATE STUFF HERE

	// Assign rate events
	localEvents.addEvent(RATE_EVENT,currentHost->wormTotalDeathRate,currentHost);
	localEvents.addEvent(RATE_EVENT,currentHost->hostInfectionRate,currentHost);

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

	return true;
}

// Is it a new worm or dead worm?
bool CRealization::wormBirthDeathResponse(Event& currentEvent)
{

	// INSERT WORM BIRTH AND DEATH STUFF HERE

	return true;
}

// Update the freeliving worm populations deterministically
bool CRealization::wormFreelivingResponse(Event& currentEvent)
{

	// INSERT WORM FREELIVING STUFF HERE

	return true;
}

// Apply chemotherapy
bool CRealization::chemotherapyResponse(Event& currentEvent)
{

	// INSERT CHEMOTHERAPY STUFF HERE

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
