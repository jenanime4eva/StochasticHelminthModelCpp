/*
 * CRealization.cpp
 *
 *  Created on: 12 Dec 2014
 *      Author: jtrusc
 */

#include "CRealization.h"
#include "CSimulator.h"

CRealization::CRealization(CSimulator* currentOwner)
{
	owner = currentOwner;
}

CRealization::~CRealization()
{
	// TODO Auto-generated destructor stub
}

// set up the realization.
bool CRealization::initialize()
{
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

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,95.0,NULL);

	return true;
}

// run the realization.
bool CRealization::run()
{
	// Get the start time.
	double time = owner->startYear;
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
			case TERMINATE:
				// nothing to do here. Perhaps shouldn't even switch on it?
				break;
			default:
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

// just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{
	owner->logStream << "\nCurrent ages at time: " << currentEvent.time << "\n";
	for(int i=0;i<nHosts;i++)
	{
		owner->logStream << currentEvent.time - hostPopulation[i]->birthDate << "\n";
	}
	owner->logStream << std::flush;
}




