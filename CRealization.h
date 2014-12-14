/*
 * CRealization.h
 *
 *  This class is encapsulate all the data involved in setting up and running a realisation.
 *  The purpose is that realizations will then be independent of each other and therefore easily parallelizable.
 */

#ifndef CREALIZATION_H_
#define CREALIZATION_H_

#include "CHost.h"
#include "CPreDetEventQueue.h"

class CRealization
{
public:
	CRealization(class CSimulator* currentOwner);
	virtual ~CRealization();

	bool initialize();
	bool run();

	bool hostDeathResponse(Event& currentEvent);
	void debugEventResponse(Event& currentEvent);

	// members...
	class CSimulator* owner;
	int nHosts;
	CHost** hostPopulation; // Array of host population.

	// realization event queue.
	CPreDetEventQueue localEvents;
};

#endif /* CREALIZATION_H_ */
