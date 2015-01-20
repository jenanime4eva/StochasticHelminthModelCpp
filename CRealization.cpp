/*
 * CRealization.cpp
 *
 *  Created on: 12 Dec 2014
 *      Author: jtrusc
 */

#include "CRealization.h"
#include "CSimulator.h"

CRealization::CRealization()
{
	hostPopulation = NULL;
	owner = NULL;
	nHosts = 0;
	surveyResultsArray = NULL;
}

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

// set up the realization.
bool CRealization::initialize(CSimulator* currentOwner)
{
	owner = currentOwner; // set the pointer to the owner.

	nHosts = owner->nHosts;

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;
		// do lifespan stuff...
		double lifespan = owner->drawLifespan();
		hostPopulation[i]->birthDate = -owner->myRand()*lifespan;
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);
	}

	// add run termination point.
	localEvents.addEvent(TERMINATE,owner->startYear+owner->nYears,NULL);

	// set up results collection for this realization.
	if(owner->surveyResultTimes!=NULL)
	{
		// there are survey results to collect.
		// for each time point set up an array of surveyResultData structures the length of the population size.
		surveyResultsArray = new surveyResultData*[owner->surveyResultTimesLength];
		for(int i=0;i<owner->surveyResultTimesLength;i++)
		{
			surveyResultsArray[i] = new surveyResultData[nHosts];
			localEvents.addEvent(SURVEY_EVENT,owner->surveyResultTimes[i],surveyResultsArray[i]); // add a pointer to the array that's going to hold the data.
		}
	}

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,95.0,NULL);

	return true;
}

// run the realization.
bool CRealization::run()
{
	// Get the start time.
	//double time = owner->startYear;
	Event currentEvent;

	do
	{
		// get next preDet Event.
		localEvents.popEvent(currentEvent);

		// do stochastic events up to the time of the next event.

		// do the event.
		switch (currentEvent.type)
		{
			case HOST_DEATH:
				// create new host here...
				hostDeathResponse(currentEvent);
				break;
			case DEBUG_EVENT:
				debugEventResponse(currentEvent);
				break;
			case SURVEY_EVENT:
				surveyResultResponse(currentEvent);
				break;
			case TERMINATE:
				// nothing to do here. Perhaps shouldn't even switch on it?
				break;
			default:
				owner->logStream << "Event number " << currentEvent.type << " not known.\n" << std::flush;
				break;
		}
	}while(currentEvent.type!=TERMINATE);

	return true;
}

//Respond to host death.
bool CRealization::hostDeathResponse(Event& currentEvent)
{
	// rejuvenate host...
	CHost* currentHost = (CHost*) currentEvent.subject;
	currentHost->birthDate = currentEvent.time;
	double lifespan = owner->drawLifespan();
	currentHost->deathDate = currentEvent.time + lifespan;

	// assign death event.
	localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

	return true;
}

bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// collect data from each host individual.
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].age = currentEvent.time - hostPopulation[i]->birthDate;
	}
	return true;
}

// just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{
}




