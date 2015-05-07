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
#include "CWorm.h"
#include "CPreDetEventQueue.h"

// Structure to contain data from individuals in survey results
struct surveyResultData
{
	double age;
	int nFemaleWorms;
	int NtotalWorms;
};

class CRealization
{
public:
	CRealization();
	virtual ~CRealization();

	bool initialize(class CSimulator* currentOwner);
	bool run();

	// Predetermined event responses
	bool hostDeathResponse(Event& currentEvent);
	bool surveyResultResponse(Event& currentEvent);
	void debugEventResponse(Event& currentEvent);

	// MEMBERS
	class CSimulator* owner;
	int nHosts;
	int femaleWorms, totalWorms;
	CHost** hostPopulation; // Array of host population
	CWorm** femaleWormNumbers; // Array of female worm numbers for each host
	CWorm** totalWormNumbers; // Array of total worm numbers for each host

	surveyResultData** surveyResultsArray;

	// Realization event queue
	CPreDetEventQueue localEvents;
};

#endif /* CREALIZATION_H_ */
