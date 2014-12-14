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

// structure to contain data from individuals in survey results.
struct surveyResultData
{
	double age;
};

class CRealization
{
public:
	CRealization();
	virtual ~CRealization();

	bool initialize(class CSimulator* currentOwner);
	bool run();

	// Pre-determined event responses...
	bool hostDeathResponse(Event& currentEvent);
	bool surveyResultResponse(Event& currentEvent);
	void debugEventResponse(Event& currentEvent);

	////////////////////////////////////////////////////////
	// members...
	class CSimulator* owner;
	int nHosts;
	CHost** hostPopulation; // Array of host population.

	//int surveyResultsArrayLength;  // got this in CSimulator.
	surveyResultData** surveyResultsArray;

	// realization event queue.
	CPreDetEventQueue localEvents;
};

#endif /* CREALIZATION_H_ */
