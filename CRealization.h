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
#include <fstream>
#include <random>
#include <string>

using namespace std;

// Structure to contain data from individuals in survey results
struct surveyResultData
{
	double age;
	double sumFemaleWorms;
	double sumInfantFemaleWorms;
	double sumPreSACFemaleWorms;
	double sumSACFemaleWorms;
	double sumAdultFemaleWorms;
	int infantNumber;
	int preSACNumber;
	int SACNumber;
	int adultNumber;
};

class CRealization
{
public:
	CRealization();
	virtual ~CRealization();

	bool initialize(class CSimulator* currentOwner,int repNo);
	bool run(int repNo);

	// Predetermined event responses
	bool hostDeathResponse(Event& currentEvent);
	bool hostChemoResponse(Event& currentEvent);
	bool surveyResultResponse(Event& currentEvent);
	void debugEventResponse(Event& currentEvent);

	// Other functions
	void calculateEventRates();
	void doEvent();
	void doFreeliving(double ts);

	// Members
	class CSimulator* owner;
	int nHosts;
	int runs;
	double nYears;
	CHost** hostPopulation; // Array of host population
	double initialWormNumber;
	double tinyIncrement;
	double freeliving;
	double* hostTotalWorms;
	double* hostFemaleWorms;
	double* productiveFemaleWorms;
	double* eggsOutputPerHost;
	double eggsProductionRate;
	double* hostInfectionRate;
	int* contactAgeGroupIndex;
	int* treatmentAgeGroupIndex;
	double* mu;
	double* rates;
	int ratesLength;
	double* nextStepCompare;
	int nextStepCompareLength;
	double* newNextStepCompare;
	int newNextStepCompareLength;
	double* chemoTimes;
	double* outTimes;
	int surveyLength;
	int treatLength;
	double ageingInterval;
	double maxStep;
	double ts;
	double timeNow;

	surveyResultData** surveyResultsArrayPerRun;

	// Realization event queue
	CPreDetEventQueue localEvents;
};

#endif /* CREALIZATION_H_ */
