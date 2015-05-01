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
	femaleWormNumbers = NULL;
	totalWormNumbers = NULL;
	owner = NULL;
	nHosts = 0;
	femaleWorms = 0;
	totalWorms = 0;
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

	if(femaleWormNumbers!=NULL)
	{
		// Female worms for each host
		for(i=0;i<nHosts;i++)
			delete femaleWormNumbers[i];

		delete[] femaleWormNumbers;
	}

	if(totalWormNumbers!=NULL)
	{
		// Female worms for each host
		for(i=0;i<nHosts;i++)
			delete totalWormNumbers[i];

		delete[] totalWormNumbers;
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
		//genunf(double low, double high) generates uniform real between low and high
		hostPopulation[i]->birthDate = -owner->myRandUni()*lifespan;
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);
	}

	// Set up female and total worm numbers array for each of the hosts
	femaleWormNumbers = new CWorm* [nHosts];
	totalWormNumbers = new CWorm* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		femaleWormNumbers[i] = new CWorm;
		totalWormNumbers[i] = new CWorm;
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
	localEvents.addEvent(DEBUG_EVENT,145.0,NULL);

	return true;
}

// Run the realisation
bool CRealization::run()
{
	Event currentEvent;

	//TODO: Add chemotherapy event

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
			//case DEBUG_EVENT:
				//debugEventResponse(currentEvent);
				//break;
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

//Respond to host death
bool CRealization::hostDeathResponse(Event& currentEvent)
{
	// Rejuvenate host
	CHost* currentHost = (CHost*) currentEvent.subject;
	currentHost->birthDate = currentEvent.time;
	double lifespan = owner->drawLifespan();
	currentHost->deathDate = currentEvent.time + lifespan;

	// Assign death event
	localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

	return true;
}

bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// Collect data from each host individual
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	// For looking at the host ages across time
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].age = currentEvent.time - hostPopulation[i]->birthDate;
	}

	/*
	// For looking at female worm burdens for individual hosts across time
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].nFemaleWorms = femaleWormNumbers[i]->femaleWorms;
	}
	*/

	return true;
}

/*
// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{
}
*/
