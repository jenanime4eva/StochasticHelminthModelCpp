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
	hostTotalWorms = new double[nHosts];
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
		double communityBurnIn = 1000.0; // 1000 years should be plenty of time to equilibrate the population
		while(hostPopulation[i]->deathDate < communityBurnIn)
		{
			double newlifespan = owner->drawLifespan();
			hostPopulation[i]->birthDate = hostPopulation[i]->deathDate;
			hostPopulation[i]->deathDate = hostPopulation[i]->deathDate + newlifespan;
		}

		// Updated birth and death dates
		hostPopulation[i]->birthDate = hostPopulation[i]->birthDate - communityBurnIn;
		hostPopulation[i]->deathDate = hostPopulation[i]->deathDate - communityBurnIn;

		// Host age
		int age = (int) floor(-hostPopulation[i]->birthDate); // 0 - birthDate (as time at this point is zero)

		// Work out hostContactIndex and hostTreatIndex
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

	// Get minimum value of the nextStepCompare array
	nextStepCompare[0] = nextOutTime;
	nextStepCompare[1] = timeNow+maxStep;
	nextStepCompare[2] = nextChemoTime;
	nextStepCompare[3] = nextAgeTime;
	double nextStep = owner->min(nextStepCompare,4);

	int outCount = 0;

	do
	{
		//printf("1st host totalworms femaleWorms: %f %f\n",hostPopulation[0]->totalWorms,hostPopulation[0]->femaleWorms);

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

			// Ageing and death
			if(timeBarrier>nextAgeTime)
			{
				nextAgeTime = nextAgeTime + ageingInterval;
			}

			// Chemotherapy
			if(timeBarrier>nextChemoTime)
			{
				chemoTimes[nextChemoIndex] = nYears + 10; // The 10 is to make sure chemo isn't been done at maxTime
				nextChemoIndex = owner->indexSmallestElement(chemoTimes,treatLength);
				nextChemoTime = chemoTimes[nextChemoIndex];
			}

			// Output
			if(timeBarrier>nextOutTime)
			{
				outCount++; // Add one to outCount
				outTimes[nextOutIndex] = nYears + 10;
				nextOutIndex = owner->indexSmallestElement(outTimes,surveyLength);
				nextOutTime = outTimes[nextOutIndex];
			}

			 // Get new minimum value of the nextStepCompare array
			timeNow = nextStep;
			//printf("%f\n",timeNow);
			newNextStepCompare[0] = nextOutTime;
			newNextStepCompare[1] = timeNow+maxStep;
			newNextStepCompare[2] = nextChemoTime;
			newNextStepCompare[3] = nextAgeTime;
			nextStep = owner->min(newNextStepCompare,4);
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
	//printf("lifespan %f\n",lifespan);

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

	// New host age
	int age = (int) floor(currentEvent.time - currentHost->birthDate);

	// Update hostContactIndex and hostTreatIndex
	currentHost->hostContactIndex = owner->contactAgeGroupIndex()[age];
	currentHost->hostTreatIndex = owner->treatmentAgeGroupIndex()[age];

	// Assign death event
	localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

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
		// These reset after each run iteration
		outputArray[rIndex].sumFemaleWorms = 0;
		outputArray[rIndex].sumInfantFemaleWorms = 0;
		outputArray[rIndex].sumPreSACFemaleWorms = 0;
		outputArray[rIndex].sumSACFemaleWorms = 0;
		outputArray[rIndex].sumAdultFemaleWorms = 0;
		outputArray[rIndex].infantNumber = 0;
		outputArray[rIndex].preSACNumber = 0;
		outputArray[rIndex].SACNumber = 0;
		outputArray[rIndex].adultNumber = 0;

		// Female worms for each host
		for(int i=0;i<nHosts;i++)
		{
			// Get current age of host
			double currentHostAge = currentEvent.time - hostPopulation[i]->birthDate;

			// Look at whole population
			outputArray[rIndex].sumFemaleWorms += hostPopulation[i]->femaleWorms;

			// Look at infants
			bool infants = (currentHostAge >= birthCutoff) && (currentHostAge < infantCutoff);

			if(infants)
			{
				outputArray[rIndex].sumInfantFemaleWorms += hostPopulation[i]->femaleWorms;
				outputArray[rIndex].infantNumber++; // Add one to infant number
			}

			// Look at pre-SAC
			bool preSAC = (currentHostAge >= infantCutoff) && (currentHostAge < preSACCutoff);

			if(preSAC)
			{
				outputArray[rIndex].sumPreSACFemaleWorms += hostPopulation[i]->femaleWorms;
				outputArray[rIndex].preSACNumber++; // Add one to pre-SAC number
			}

			// Look at SAC
			bool SAC = (currentHostAge >= preSACCutoff) && (currentHostAge < SACCutoff);

			if(SAC)
			{
				outputArray[rIndex].sumSACFemaleWorms += hostPopulation[i]->femaleWorms;
				outputArray[rIndex].SACNumber++; // Add one to SAC number
			}

			// Look at adults
			bool adults = (currentHostAge >= SACCutoff) && (currentHostAge < adultCutoff);

			if(adults)
			{
				outputArray[rIndex].sumAdultFemaleWorms += hostPopulation[i]->femaleWorms;
				outputArray[rIndex].adultNumber++; // Add one to adult number
			}
		}

		// Demography test
		printf("infant, pre-SAC, SAC and adult count is %d %d %d %d\n",outputArray[rIndex].infantNumber,outputArray[rIndex].preSACNumber,outputArray[rIndex].SACNumber,outputArray[rIndex].adultNumber);
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
		// Store individual host total worms into the hostTotalWorms array for use in the multiNomBasic2 function
		for (int i=0;i<nHosts;i++)
		{
			hostTotalWorms[i] = hostPopulation[i]->totalWorms;
		}

		// Dead worm
		int deathIndex = owner->multiNomBasic2(hostTotalWorms,nHosts,owner->myRandUniform());

		double femaleMaleWormProportion = (hostPopulation[deathIndex]->femaleWorms) / (hostPopulation[deathIndex]->totalWorms);

		bool femaleWormCondition = owner->myRandUniform() < femaleMaleWormProportion;

		// Is this worm female?
		if(femaleWormCondition)
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

		bool newFemaleWormCondition = owner->myRandUniform() < 0.5;

		// New female worm
		if(newFemaleWormCondition)
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
		bool noMales = hostPopulation[i]->totalWorms == hostPopulation[i]->femaleWorms;

		if (noMales)
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
