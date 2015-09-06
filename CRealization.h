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
	double hostAge;
	int femaleWorms;
	int infantFemaleWorms;
	int preSACFemaleWorms;
	int SACFemaleWorms;
	int adultFemaleWorms;
	int infantNumber;
	int preSACNumber;
	int SACNumber;
	int adultNumber;
	int totalHostNumber;
};

class CRealization
{
public:
	CRealization();
	virtual ~CRealization();

	bool initialize(class CSimulator* currentOwner,int repNo);
	bool run(int repNo);

	// Predetermined event responses
	bool surveyResultResponse(Event& currentEvent);

	// Other functions
	void calculateEventRates(int sumTotalWorms);
	void doEvent(double* hostTotalWorms);
	double doFreeliving(double ts,int freeliving);
	void doDeath(double timeNow);
	void doChemo();

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
	int* productiveFemaleWorms;
	double* eggsOutputPerHost;
	double eggsProductionRate;
	double* hostInfectionRate;
	int* contactAgeGroupIndex;
	int* treatmentAgeGroupIndex;
	double* mu;
	double* rates;
	double* compareArray;
	double* chemoTimes;
	double* outTimes;
	int ratesLength;
	int surveyLength;
	int treatLength;
	double ageingInterval;
	double maxStep;
	double ts;
	double timeNow;
	double birthCutoff;
	double infantCutoff;
	double preSACCutoff;
	double SACCutoff;
	double adultCutoff;

	int counter; // Global variable counter

	surveyResultData** surveyResultsArrayPerRun;

	// Realization event queue
	CPreDetEventQueue localEvents;
};

#endif /* CREALIZATION_H_ */
