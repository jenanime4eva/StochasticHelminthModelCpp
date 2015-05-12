/*
 * File Name: CRealization.h
 *
 * Created on: 14 April 2015
 * Author: Jie Yang
 *
 * This class is to encapsulate all the data involved in setting up and running a realisation
 * The purpose is that realisations will then be independent of each other and therefore easily parallelisable
 *
 */

#ifndef CREALIZATION_H_
#define CREALIZATION_H_

#include "CHost.h"
#include "CPreDetEventQueue.h"

// Structure to contain data from individuals in survey results
struct surveyResultData
{
	double age;
	double nFemaleWorms;
	double nTotalWorms;
};

class CRealization
{
public:
	CRealization();
	virtual ~CRealization();

	bool initialize(class CSimulator* currentOwner);
	bool run();

	// Predetermined event responses
	bool calculateEventRatesResponse(Event& currentEvent);
	bool hostDeathResponse(Event& currentEvent);
	bool wormBirthDeathResponse(Event& currentEvent);
	bool wormFreelivingResponse(Event& currentEvent);
	bool chemotherapyResponse(Event& currentEvent);
	bool surveyResultResponse(Event& currentEvent);
	void debugEventResponse(Event& currentEvent);

	// MEMBERS
	class CSimulator* owner;
	int nHosts;
	CHost** hostPopulation; // Array of host population

	surveyResultData** surveyResultsArray;

	// Realization event queue
	CPreDetEventQueue localEvents;
};

#endif /* CREALIZATION_H_ */
