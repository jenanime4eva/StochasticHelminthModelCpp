/*
 * File Name: CRealization.cpp
 *
 *  Created on: 14 April 2015
 *  Author: Jie Yang
 *
 */

#include "CRealization.h"
#include "CSimulator.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <vector.h>
#include "randlib.h"
#include <time.h>

// Class constructor
CRealization::CRealization()
{
	hostPopulation = NULL;
	owner = NULL;
	nHosts = 0;
	runs = 0;
	nYears = 0;
	initialWormNumber = 11.5; // Choose this number such that initial conditions are at the equilibrium
	tinyIncrement = 0.01;
	freeliving = 0;
	surveyResultsArrayPerRun = NULL;
	hostTotalWorms = NULL;
	hostFemaleWorms = NULL;
	productiveFemaleWorms = NULL;
	eggsOutputPerHost = NULL;
	eggsProductionRate = 0;
	hostInfectionRate = NULL;
	contactAgeGroupIndex = NULL;
	treatmentAgeGroupIndex = NULL;
	mu = NULL;
	rates = NULL;
	ratesLength = 0;
	nextStepCompare = NULL;
	nextStepCompareLength = 4;
	newNextStepCompare = NULL;
	newNextStepCompareLength = 4;
	chemoTimes = NULL;
	outTimes = NULL;
	surveyLength = 0;
	treatLength = 0;
	ageingInterval = 0.25; // Checks ages every quarter year
	maxStep = (double) 1/52; // Time steps never exceed this for deterministic update of freeliving populations
	ts = 0;
	timeNow = 0;
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

	if (hostTotalWorms != NULL)
		delete[] hostTotalWorms;

	if (hostFemaleWorms != NULL)
		delete[] hostFemaleWorms;

	if (productiveFemaleWorms != NULL)
		delete[] productiveFemaleWorms;

	if (eggsOutputPerHost != NULL)
		delete[] eggsOutputPerHost;

	if (hostInfectionRate != NULL)
		delete[] hostInfectionRate;

	if (contactAgeGroupIndex != NULL)
		delete[] contactAgeGroupIndex;

	if (treatmentAgeGroupIndex != NULL)
		delete[] treatmentAgeGroupIndex;

	if(mu!=NULL)
		delete[] mu;

	if (rates != NULL)
		delete[] rates;

	if (nextStepCompare != NULL)
		delete[] nextStepCompare;

	if (newNextStepCompare != NULL)
		delete[] newNextStepCompare;

	if (chemoTimes != NULL)
		delete[] chemoTimes;

	if (outTimes != NULL)
		delete[] outTimes;

	if(surveyResultsArrayPerRun!=NULL)
	{
		for(i=0;i<owner->surveyResultTimesLength;i++)
			delete[] surveyResultsArrayPerRun[i];

		delete[] surveyResultsArrayPerRun;
	}
}

// Set up the realisation
bool CRealization::initialize(CSimulator* currentOwner,int repNo)
{
	// Set the pointer to the owner
	owner = currentOwner;

	nHosts = owner->nHosts; // Number of hosts
	runs = owner->nRepetitions; // Number of runs
	nYears = owner->nYears; // Number of years to run
	surveyLength = owner->surveyResultTimesLength; // Length of survey times array
	treatLength = owner->treatmentTimesLength; // Length of treatment times array

	// Set up some arrays required later
	// Total and female host worm arrays
	hostTotalWorms = new double[nHosts];
	hostFemaleWorms = new double[nHosts];
	productiveFemaleWorms = new double[nHosts];
	eggsOutputPerHost = new double[nHosts];
	hostInfectionRate = new double[nHosts];

	// Set up some arrays for use later
	mu = new double[nHosts];
	ratesLength = nHosts+1;
	rates = new double[ratesLength];
	nextStepCompare = new double[nextStepCompareLength];
	newNextStepCompare = new double[newNextStepCompareLength];
	chemoTimes = new double[treatLength];
	outTimes = new double[surveyLength];

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Initialise total and female worms per individual host
		// Force of infection
		hostPopulation[i]->si = owner->myRandGamma(owner->k,owner->k); // location = 1/scale (here scale = 1/k, so location = k); shape parameter is the parameter k
		mu[i] = 2.0*initialWormNumber*hostPopulation[i]->si; // These values chosen as initial conditions (initial conditions should be the equilibrium)
		hostPopulation[i]->totalWorms = owner->myRandPoisson(mu[i]); // Total worms
		hostPopulation[i]->femaleWorms = owner->myRandBinomial(hostPopulation[i]->totalWorms,0.5); // Half of worm population should be female

		// Initial freeliving worms, don't know what this number should be
		freeliving = 4.0;

		// Distribute birth dates such that at time zero have the right distribution of ages
		// Give a sample of ages at the start of the simulation
		// For an exponential distribution, the survival function S(t) = exp(-t/mu)
		// So to sample from this distribution, generate a random number on 0 1 and then invert this function

		// Host life span
		double lifespan = owner->drawLifespan();

		// Host birth and death date
		hostPopulation[i]->birthDate = -owner->myRandUniform()*lifespan;
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;

		// Equilibrate the population first, in the absence of understanding how to generate it in the first place
		double communityBurnIn = 1000.0;
		while(hostPopulation[i]->deathDate < communityBurnIn)
		{
			double newlifespan = owner->drawLifespan();
			hostPopulation[i]->birthDate = hostPopulation[i]->deathDate;
			hostPopulation[i]->deathDate = hostPopulation[i]->deathDate + newlifespan;
		}

		// Updated birth and death dates
		hostPopulation[i]->birthDate = hostPopulation[i]->birthDate - communityBurnIn;
		hostPopulation[i]->deathDate = hostPopulation[i]->deathDate - communityBurnIn;

		// Work out hostContactIndex and hostTreatIndex
		hostPopulation[i]->hostAge = -hostPopulation[i]->birthDate; // 0 - birthDate (as time at this point is zero)
		int age = (int) floor(hostPopulation[i]->hostAge);
		hostPopulation[i]->hostContactIndex = owner->contactAgeGroupIndex()[age];
		hostPopulation[i]->hostTreatIndex = owner->treatmentAgeGroupIndex()[age];

		// Add host death events
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);
	}

	// Add treatment events
	for(int j=0;j<treatLength;j++)
	{
		localEvents.addEvent(CHEMOTHERAPY,owner->treatmentTimes[j],NULL);
	}

	// Set up results collection for this realisation
	if(owner->surveyResultTimes!=NULL)
	{
		// There are survey results to collect

		// For each time point set up an array of surveyResultData structures the length of the number of realisations
		surveyResultsArrayPerRun = new surveyResultData*[surveyLength];
		for(int i = 0; i < surveyLength; i++)
		{
			// For each time point allocate nRepetitions worth of space
			surveyResultsArrayPerRun[i] = new surveyResultData[runs];
			// Add a pointer to the arrays that are going to hold the data
			localEvents.addEvent(SURVEY_EVENT,owner->surveyResultTimes[i],surveyResultsArrayPerRun[i]);
		}
	}

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,nYears - 5.0,NULL);

	return true;
}

// Run the realisation
bool CRealization::run(int repNo)
{
	//printf("Run number: %d\n",repNo); // Test flag

	Event currentEvent;

	// Some initial time-related variables
	timeNow = 0;
	double FLlast = timeNow;

	for(int i=0;i<surveyLength;i++)
	{
		outTimes[i] = owner->surveyResultTimes[i];
	}

	int nextOutIndex = 0;
	double nextOutTime = outTimes[nextOutIndex]; // First survey output time
	double nextAgeTime = ageingInterval;

	for(int j=0;j<treatLength;j++)
	{
		chemoTimes[j] = owner->treatmentTimes[j];
	}

	int nextChemoIndex = 0;
	double nextChemoTime = chemoTimes[nextChemoIndex]; // First chemotherapy time

	nextStepCompare[0] = nextOutTime;
	nextStepCompare[1] = timeNow+maxStep;
	nextStepCompare[2] = nextChemoTime;
	nextStepCompare[3] = nextAgeTime;
	double nextStep = owner->min(nextStepCompare,4); // Get minimum value of the nextStepCompare array

	int outCount = 0;

	do
	{
		// Calculate the event rates (host infection rate and worm total death rate)
		calculateEventRates();

		// Calculate sum of rates
		double sumRates = owner->sumArray(rates,ratesLength);

		// Calculate the time step
		double tstep = owner->myRandExponential(sumRates);

		if( (timeNow+tstep) < nextStep )
		{
			timeNow = timeNow + tstep;

			// Enact an event
			doEvent();
		}
		else // (timeNow+tstep) >= nextStep
		{
			// Time step
			ts = nextStep - FLlast;

			// Update the freeliving worm populations (worms found in the environment) deterministically
			doFreeliving(ts);

			FLlast = nextStep;

			// Get next predetermined event
			localEvents.popEvent(currentEvent);

			// Do the event
			switch (currentEvent.type)
			{
				// Death of host
				case HOST_DEATH:
					hostDeathResponse(currentEvent);
					break;

				// Apply treatment
				case CHEMOTHERAPY:
					hostChemoResponse(currentEvent);
					break;

				// Any events to debug?
				case DEBUG_EVENT:
					debugEventResponse(currentEvent);
					break;

				// What to output from the simulation
				case SURVEY_EVENT:
					surveyResultResponse(currentEvent);
					break;

				default:
					owner->logStream << "Event number " << currentEvent.type << " not known.\n"
					<< "Look at #define values in CPreDetEventQueue.h for event definition.\n" << std::flush;
					break;
			}

			double timeBarrier = nextStep + 0.001;

			if(timeBarrier>nextAgeTime)
			{
				nextAgeTime = nextAgeTime + ageingInterval;
			}

			if(timeBarrier>nextChemoTime)
			{
				chemoTimes[nextChemoIndex] = nYears + 10; // The 10 is to make sure chemo isn't been done at maxTime
				nextChemoIndex = owner->indexSmallestElement(chemoTimes,treatLength);
				nextChemoTime = chemoTimes[nextChemoIndex];
			}

			if(timeBarrier>nextOutTime)
			{
				outCount++;
				outTimes[nextOutIndex] = nYears + 10;
				nextOutIndex = owner->indexSmallestElement(outTimes,surveyLength);
				nextOutTime = outTimes[nextOutIndex];
			}

			timeNow = nextStep;
			newNextStepCompare[0] = nextOutTime;
			newNextStepCompare[1] =timeNow+maxStep;
			newNextStepCompare[2] = nextChemoTime;
			newNextStepCompare[3] = nextAgeTime;
			nextStep = owner->min(newNextStepCompare,4); // Get new minimum value of the nextStepCompare array

		} // End of predetermined event block
	} while((timeNow<nYears) && (outCount<surveyLength)); // Do events until both nYears and the last survey year is reached

	return true;
}


//////////////////////////////////////////////////////////////////////////////////
/// PREDETERMINED EVENT RESPONSES

// Respond to host death
bool CRealization::hostDeathResponse(Event& currentEvent)
{
	// Rejuvenate host
	CHost* currentHost = (CHost*) currentEvent.subject;

	// Host life span
	double lifespan = owner->drawLifespan();

	// Set birth date to now
	// Put birth slightly in the past to ensure age is just positive for catergorisation
	currentHost->birthDate = currentEvent.time - 0.001;

	// Calculate new death dates
	currentHost->deathDate = currentEvent.time + lifespan;

	// Calculate new force of infection (FOI)
	currentHost->si = owner->myRandGamma(owner->k,owner->k);

	// Kill off all their worms
	currentHost->totalWorms = 0.0;
	currentHost->femaleWorms = 0.0;

	// Update hostContactIndex and hostTreatIndex
	currentHost->hostAge = currentEvent.time - currentHost->birthDate;
	int age = (int) floor(currentHost->hostAge);
	currentHost->hostContactIndex = owner->contactAgeGroupIndex()[age];
	currentHost->hostTreatIndex = owner->treatmentAgeGroupIndex()[age];

	// Assign death event
	//localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

	return true;
}

// Chemotherapy
bool CRealization::hostChemoResponse(Event& currentEvent)
{
	for(int i=0;i<nHosts;i++)
	{
		// Coverage level for each host
		double individualCoverage = owner->coverage[hostPopulation[i]->hostTreatIndex];

		// Condition for those treated, randomly chosen
		bool individualTreated = owner->myRandUniform() < individualCoverage;

		// If individualTreated condition is TRUE
		if (individualTreated)
		{
			// How many worms to die?
			double totalW = hostPopulation[i]->totalWorms;
			double femaleW = hostPopulation[i]->femaleWorms;
			double maleW = totalW - femaleW;

			double drugEfficacy = owner->drugEff;

			double maleToDie = owner->myRandBinomial(maleW,drugEfficacy);
			double femaleToDie = owner->myRandBinomial(femaleW,drugEfficacy);

			hostPopulation[i]->totalWorms = totalW - maleToDie - femaleToDie;
			hostPopulation[i]->femaleWorms = femaleW - femaleToDie;
		}
	}

	return true;
}

// What to output from the simulation
bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// Collect data from each run
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	// Treatment age cutoffs
	double birthCutoff = owner->treatmentBreaks[0];
	double infantCutoff = owner->treatmentBreaks[1];
	double preSACCutoff = owner->treatmentBreaks[2];
	double SACCutoff = owner->treatmentBreaks[3];
	double adultCutoff = owner->treatmentBreaks[4];

	// Collect data from each run
	for(int rIndex=0;rIndex<runs;rIndex++)
	{
		// Set up some variables that resets after each iteration
		double sumFemaleWormsPerRun = 0.0;
		double sumInfantFemaleWormsPerRun = 0.0;
		double sumPreSACFemaleWormsPerRun = 0.0;
		double sumSACFemaleWormsPerRun = 0.0;
		double sumAdultFemaleWormsPerRun = 0.0;
		int infantCount = 0; // Infant number counter
		int preSACCount = 0; // Pre-SAC number counter
		int SACCount = 0; // SAC number counter
		int adultCount = 0; // Adult number counter

		// Female worms for each host
		for(int i=0;i<nHosts;i++)
		{
			// Get current age of host
			double currentHostAge = hostPopulation[i]->hostAge;

			// Looks at whole population
			sumFemaleWormsPerRun += hostPopulation[i]->femaleWorms;

			// Look at infants
			if((currentHostAge >= birthCutoff) && (currentHostAge < infantCutoff) )
			{
				sumInfantFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
				infantCount++; // Add one to infant number counter
			}

			// Look at pre-SAC
			if((currentHostAge >= infantCutoff) && (currentHostAge < preSACCutoff) )
			{
				sumPreSACFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
				preSACCount++; // Add one to pre-SAC number counter
			}

			// Look at SAC
			if((currentHostAge >= preSACCutoff) && (currentHostAge < SACCutoff) )
			{
				sumSACFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
				SACCount++; // Add one to SAC number counter
			}

			// Look at adults
			if((currentHostAge >= SACCutoff) && (currentHostAge < adultCutoff) )
			{
				sumAdultFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
				adultCount++; // Add one to adult number counter
			}
		}

		//printf("infant, pre-SAC, SAC and adult count is %d %d %d %d\n",infantCount,preSACCount,SACCount,adultCount);
		//printf("Female numbers %f %f %f %f\n",sumInfantFemaleWormsPerRun,sumPreSACFemaleWormsPerRun,sumSACFemaleWormsPerRun,sumAdultFemaleWormsPerRun);

		// Whole population
		outputArray[rIndex].sumFemaleWorms = sumFemaleWormsPerRun;

		// Infants
		outputArray[rIndex].sumInfantFemaleWorms = sumInfantFemaleWormsPerRun;
		outputArray[rIndex].infantNumber = infantCount;

		// Pre-SAC
		outputArray[rIndex].sumPreSACFemaleWorms = sumPreSACFemaleWormsPerRun;
		outputArray[rIndex].preSACNumber = preSACCount;

		// SAC
		outputArray[rIndex].sumSACFemaleWorms = sumSACFemaleWormsPerRun;
		outputArray[rIndex].SACNumber = SACCount;

		// Adults
		outputArray[rIndex].sumAdultFemaleWorms = sumAdultFemaleWormsPerRun;
		outputArray[rIndex].adultNumber = adultCount;
	}

	return true;
}

// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{

}


//////////////////////////////////////////////////////////////////////////////////
/// OTHER FUNCTIONS

// Calculate the event rates
void CRealization::calculateEventRates()
{
	// Reset these variables to zero after each iteration
	double sumTotalWorms = 0.0;

	for (int i=0;i<nHosts;i++)
	{
		hostInfectionRate[i] = freeliving*hostPopulation[i]->si*owner->betaValues[hostPopulation[i]->hostContactIndex];
		rates[i] = hostInfectionRate[i];
		sumTotalWorms += hostPopulation[i]->totalWorms; // Sum of total worms for all hosts
	}
	double wormTotalDeathRate = owner->sigma*sumTotalWorms;
	rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array
}

// Enact an event
void CRealization::doEvent()
{
	// Determine if dead worm or new worm
	int event = owner->multiNomBasic2(rates,ratesLength,owner->myRandUniform());

	if(event==ratesLength-1)
	{
		for (int i=0;i<nHosts;i++)
		{
			hostTotalWorms[i] = hostPopulation[i]->totalWorms;
			hostFemaleWorms[i] = hostPopulation[i]->femaleWorms;
		}

		// Dead worm
		int deathIndex = owner->multiNomBasic2(hostTotalWorms,nHosts,owner->myRandUniform());

		// Is this worm female?
		if(owner->myRandUniform() < (owner->myRandUniform(),hostPopulation[deathIndex]->femaleWorms/hostPopulation[deathIndex]->totalWorms))
		{
			// Remove a worm from female worms
			hostPopulation[deathIndex]->femaleWorms = hostPopulation[deathIndex]->femaleWorms - 1.0;
		}

		// Remove worm from total worms
		hostPopulation[deathIndex]->totalWorms = hostPopulation[deathIndex]->totalWorms - 1.0;
	}
	else
	{
		// New worm in total worms
		hostPopulation[event]->totalWorms = hostPopulation[event]->totalWorms + 1.0;

		// New female worm
		if(owner->myRandUniform() < 0.5)
		{
			hostPopulation[event]->femaleWorms = hostPopulation[event]->femaleWorms + 1.0;
		}
	}

}

// Update the freeliving populations deterministically
void CRealization::doFreeliving(double ts)
{
	double sumEggsOutputPerHostRho = 0.0; // Reset this value before next iteration

	for(int i=0;i<nHosts;i++)
	{
		productiveFemaleWorms[i] = hostPopulation[i]->femaleWorms;

		// Female worms produce fertilised eggs only if there is a male worm around
		if (hostPopulation[i]->totalWorms == hostPopulation[i]->femaleWorms)
		{
			productiveFemaleWorms[i] = 0.0;
		}

		eggsOutputPerHost[i] = owner->lambda*productiveFemaleWorms[i]*exp(-productiveFemaleWorms[i]*owner->gamma);

		sumEggsOutputPerHostRho += eggsOutputPerHost[i]*owner->rhoValues[hostPopulation[i]->hostContactIndex];
	}

	eggsProductionRate = owner->psi*sumEggsOutputPerHostRho/nHosts;

	// dL/dt = K-mu*L has solution: L(0)exp(-mu*t)+K*(1-exp(-mu*t))/mu, this is exact if rate of egg production is constant in the timestep
	// Here L = freeliving, K = eggsProductionRate, mu = ReservoirDecayRate
	double expo = exp(-owner->ReservoirDecayRate*ts);
	freeliving = freeliving*expo + eggsProductionRate*(1.0-expo)/(owner->ReservoirDecayRate);

}

