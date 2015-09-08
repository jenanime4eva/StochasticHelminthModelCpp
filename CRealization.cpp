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
	compareArray = NULL;
	chemoTimes = NULL;
	outTimes = NULL;
	surveyLength = 0;
	treatLength = 0;
	ageingInterval = 0.25; // Checks ages every quarter year
	maxStep = (double) 1/52; // Time steps never exceed this for deterministic update of freeliving populations
	ts = 0;
	birthCutoff = 0;
	infantCutoff = 0;
	preSACCutoff = 0;
	SACCutoff = 0;
	adultCutoff = 0;
	counter1 = counter2 = 0; // Global variable counters for testing purposes
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

	if (compareArray != NULL)
		delete[] compareArray;

	if (chemoTimes != NULL)
		delete[] chemoTimes;

	if (outTimes != NULL)
		delete[] outTimes;

	if(surveyResultsArrayPerRun!=NULL)
	{
		for(i=0;i<surveyLength;i++)
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
	compareArray = new double[4];
	chemoTimes = new double[treatLength];
	outTimes = new double[surveyLength];

	// Treatment age cutoffs
	birthCutoff = owner->treatmentBreaks[0];
	infantCutoff = owner->treatmentBreaks[1];
	preSACCutoff = owner->treatmentBreaks[2];
	SACCutoff = owner->treatmentBreaks[3];
	adultCutoff = owner->treatmentBreaks[4];

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Initialise total and female worms per individual host
		// Force of infection
		double k = owner->k;
		hostPopulation[i]->si = owner->myRandGamma(k,k); // location = 1/scale (here scale = 1/k, so location = k); shape parameter is the parameter k
		double initialWormNumber = 11.5; // Choose this number such that initial conditions are at the equilibrium
		mu[i] = 2.0*initialWormNumber*hostPopulation[i]->si; // These values chosen as initial conditions (initial conditions should be the equilibrium)
		hostPopulation[i]->totalWorms = owner->myRandPoisson(mu[i]); // Total worms
		hostPopulation[i]->femaleWorms = owner->myRandBinomial(hostPopulation[i]->totalWorms,0.5); // Half of worm population should be female

		// Distribute birth dates such that at time zero have the right distribution of ages
		// Give a sample of ages at the start of the simulation
		// For an exponential distribution, the survival function S(t) = exp(-t/mu)
		// So to sample from this distribution, generate a random number on 0 1 and then invert this function

		// Host life span
		double lifespan = owner->drawLifespan();

		// Host birth and death date
		double randNum = owner->myRandUniform();
		hostPopulation[i]->birthDate = -randNum*lifespan;
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
		hostPopulation[i]->hostAge = -hostPopulation[i]->birthDate; // 0 - birthDate (as time at this point is zero)

		int age = (int) floor(hostPopulation[i]->hostAge);

		// Work out hostContactIndex and hostTreatIndex
		hostPopulation[i]->hostContactIndex = owner->contactAgeGroupIndex()[age];
		hostPopulation[i]->hostTreatIndex = owner->treatmentAgeGroupIndex()[age];

		// Add host death events
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);
	}

	freeliving = 4.0; // Initial freeliving worms, don't know what this number should be

	for(int j=0;j<treatLength;j++)
	{
		chemoTimes[j] = owner->treatmentTimes[j];
	}

	// Add treatment events
	for(int j=0;j<treatLength;j++)
	{
		localEvents.addEvent(CHEMOTHERAPY,chemoTimes[j],NULL);
	}

	for(int i=0;i<surveyLength;i++)
	{
		outTimes[i] = owner->surveyResultTimes[i];
	}

	// Set up results collection for this realization.
	if(outTimes!=NULL)
	{
		// There are survey results to collect
		// For each time point set up an array of surveyResultData structures the length of the population size
		surveyResultsArrayPerRun = new surveyResultData*[surveyLength];
		for(int i=0;i<surveyLength;i++)
		{
			surveyResultsArrayPerRun[i] = new surveyResultData[nHosts];
			localEvents.addEvent(SURVEY_EVENT,outTimes[i],surveyResultsArrayPerRun[i]); // Add a pointer to the array that's going to hold the data
		}
	}

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,nYears - 5.0,NULL);

	return true;
}

// Run the realisation
bool CRealization::run(int repNo)
{
	// Some initial time-related variables
	double timeNow = 0;
	double FLlast = timeNow;

	int nextOutIndex = 0;
	double nextOutTime = outTimes[nextOutIndex]; // First survey output time
	double nextAgeTime = ageingInterval;

	int nextChemoIndex = 0;
	double nextChemoTime = chemoTimes[nextChemoIndex]; // First chemotherapy time

	// Get minimum value of compareArray
	compareArray[0] = nextOutTime;
	compareArray[1] = timeNow + maxStep;
	compareArray[2] = nextChemoTime;
	compareArray[3] = nextAgeTime;

	double nextStep = owner->min(compareArray,4);

	Event currentEvent;

	do
	{
		// Calculate the event rates (host infection rate and worm total death rate)
		calculateEventRates();

		// rates is a cumulative array, what is the final entry?
		double lastRate = rates[ratesLength-1];
		//printf("lastRate %f\n",lastRate);

		// Calculate the time step
		double tstep = owner->myRandExponential(lastRate);
		//printf("tstep: %f\n",tstep);

		if( (timeNow+tstep) < nextStep )
		{
			counter1++;
			timeNow = timeNow + tstep;

			// Enact an event
			doEvent();
		}
		else // (timeNow+tstep) >= nextStep
		{
			counter2++;
			// Time step
			ts = nextStep - FLlast;
			//printf("ts: %f\n",ts);

			// Update the freeliving worm populations (worms found in the environment) deterministically
			freeliving = doFreeliving(ts,freeliving);
			//printf("freeliving %f\n",freeliving);

			FLlast = nextStep;

			// Get next predetermined event
			localEvents.popEvent(currentEvent);

			// Do the event
			switch (currentEvent.type)
			{
				// Death of host
				case HOST_DEATH:
					//printf("hostdeath currentEvent.time: %f\t",currentEvent.time);
					hostDeathResponse(currentEvent);
					nextAgeTime = nextAgeTime + ageingInterval;
					break;

				// Apply treatment
				case CHEMOTHERAPY:
					//printf("chemo currentEvent.time: %f\t",currentEvent.time);
					hostChemoResponse(currentEvent);
					chemoTimes[nextChemoIndex] = nYears + 10; // The 10 is to make sure chemo isn't been done at maxTime
					nextChemoIndex = owner->indexSmallestElement(chemoTimes,treatLength);
					nextChemoTime = chemoTimes[nextChemoIndex];
					break;

				// Any events to debug?
				case DEBUG_EVENT:
					debugEventResponse(currentEvent);
					break;

				// What to output from the simulation
				case SURVEY_EVENT:
					//printf("survey currentEvent.time: %f\t",currentEvent.time);
					surveyResultResponse(currentEvent);
					outTimes[nextOutIndex] = nYears + 10;
					nextOutIndex = owner->indexSmallestElement(outTimes,surveyLength);
					nextOutTime = outTimes[nextOutIndex];
					break;

				default:
					owner->logStream << "Event number " << currentEvent.type << " not known.\n"
					<< "Look at #define values in CPreDetEventQueue.h for event definition.\n" << std::flush;
					break;
			}

			timeNow = nextStep;

			// Get new minimum of compareArray
			compareArray[0] = nextOutTime;
			compareArray[1] = currentEvent.time + maxStep;
			compareArray[2] = nextChemoTime;
			compareArray[3] = nextAgeTime;

			nextStep = owner->min(compareArray,4);

			// TODO
			//printf("timeNow nextStep: %f %f\n",timeNow,nextStep);
			//printf("hostPopulation[0] %d\n",hostPopulation[0]->femaleWorms);

		} // End of predetermined event block
	} while(timeNow < nYears); // Do events until nYears is reached
	//printf("if %d times\n",counter1);
	//printf("else %d times",counter2);
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
	double k = owner->k;
	currentHost->si = owner->myRandGamma(k,k);

	// Kill off all their worms
	currentHost->totalWorms = 0;
	currentHost->femaleWorms = 0;

	// Update hostContactIndex and hostTreatIndex
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
		double randNum = owner->myRandUniform();
		bool individualTreated = randNum < individualCoverage;

		// If individualTreated condition is TRUE
		if (individualTreated)
		{
			// How many worms to die?
			double totalW = (double) hostPopulation[i]->totalWorms;
			double femaleW = (double) hostPopulation[i]->femaleWorms;
			double maleW = totalW - femaleW;

			double drugEfficacy = owner->drugEff;

			int maleToDie = owner->myRandBinomial(maleW,drugEfficacy);
			int femaleToDie = owner->myRandBinomial(femaleW,drugEfficacy);

			hostPopulation[i]->totalWorms = hostPopulation[i]->totalWorms - maleToDie - femaleToDie;
			hostPopulation[i]->femaleWorms = hostPopulation[i]->femaleWorms - femaleToDie;
		}
	}

	return true;
}

// What to output from the simulation
bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// Collect data from each run
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	// Loop through the hosts...
	for(int i=0;i<nHosts;i++)
	{
		// Get current age of host
		outputArray[i].hostAge = currentEvent.time - hostPopulation[i]->birthDate;
		//outputArray[i].hostAge = hostPopulation[i]->hostAge;
		double currentHostAge = outputArray[i].hostAge;

		// Look at whole population
		outputArray[i].femaleWorms = hostPopulation[i]->femaleWorms; // Number of female worms for ith host
		outputArray[i].totalHostNumber = 1;

		// Look at infants
		bool infants = (currentHostAge >= birthCutoff) && (currentHostAge < infantCutoff);

		if(infants)
		{
			outputArray[i].infantFemaleWorms = hostPopulation[i]->femaleWorms;
			outputArray[i].infantNumber = 1;
		}

		if(!infants)
		{
			outputArray[i].infantFemaleWorms = 0;
			outputArray[i].infantNumber = 0;
		}

		// Look at pre-SAC
		bool preSAC = (currentHostAge >= infantCutoff) && (currentHostAge < preSACCutoff);

		if(preSAC)
		{
			outputArray[i].preSACFemaleWorms = hostPopulation[i]->femaleWorms;
			outputArray[i].preSACNumber = 1;
		}

		if(!preSAC)
		{
			outputArray[i].preSACFemaleWorms = 0;
			outputArray[i].preSACNumber = 0;
		}

		// Look at SAC
		bool SAC = (currentHostAge >= preSACCutoff) && (currentHostAge < SACCutoff);

		if(SAC)
		{
			outputArray[i].SACFemaleWorms = hostPopulation[i]->femaleWorms;
			outputArray[i].SACNumber = 1;
		}

		if(!SAC)
		{
			outputArray[i].SACFemaleWorms = 0;
			outputArray[i].SACNumber = 0;
		}

		// Look at adults
		bool adults = (currentHostAge >= SACCutoff) && (currentHostAge < adultCutoff);

		if(adults)
		{
			outputArray[i].adultFemaleWorms = hostPopulation[i]->femaleWorms;
			outputArray[i].adultNumber = 1;
		}

		if(!adults)
		{
			outputArray[i].adultFemaleWorms = 0;
			outputArray[i].adultNumber = 0;
		}
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
	// Store individual host total worms into the hostTotalWorms array
	for (int i=0;i<nHosts;i++)
	{
		hostTotalWorms[i] = (double) hostPopulation[i]->totalWorms;
	}

	// Work out sum of the total worms array
	int sumTotalWorms = (int) owner->sumArray(hostTotalWorms,nHosts);

	// rates array has nHost entries
	// Entries 0 to nHost-1 are cumulative hostInfectionRate values, entry nHost is wormTotalDeathRate + cumulative value
	double cumul = 0; // A cumulative value
	for (int i=0;i<nHosts;i++)
	{
		double betaValue = owner->betaValues[hostPopulation[i]->hostContactIndex];
		hostInfectionRate[i] = freeliving*hostPopulation[i]->si*betaValue;
		rates[i] = cumul + hostInfectionRate[i];
		cumul = rates[i];
	}
	double sigma = owner->sigma;
	double wormTotalDeathRate = (double) (sigma*sumTotalWorms);
	rates[nHosts] = cumul + wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array
}

// Enact an event
void CRealization::doEvent()
{
	// Determine if dead worm or new worm
	// If event is equal to 0,1,2,...,nHosts-1 it's a NEW WORM, otherwise (if event is equal to nHost) it's a DEAD WORM
	double currentRand = owner->myRandUniform();
	int event = owner->multiNomBasic(rates,ratesLength,currentRand); // TODO: Check if this is working properly
	//printf("event %d\n",event);

	if(event==nHosts) // DEAD WORM
	{
		// Now let's make hostTotalWorms a cumulative array (so that multiNomBasic may be used in a bit)
		double cumulTotalWorms = 0;
		for (int i=0;i<nHosts;i++)
		{
			double TW = (double) hostPopulation[i]->totalWorms;
			hostTotalWorms[i] = cumulTotalWorms + TW;
			cumulTotalWorms = hostTotalWorms[i];
		}

		// Dead worm
		double currentRand = owner->myRandUniform();
		int deathIndex = owner->multiNomBasic(hostTotalWorms,nHosts,currentRand);
		//printf("deathIndex %d\n",deathIndex);

		// Does this host have any worms to begin with? You can't remove a worm if there weren't any to begin with!
		bool hostHasWorms = hostPopulation[deathIndex]->totalWorms != 0;

		if(hostHasWorms)
		{
			double FW = (double) hostPopulation[deathIndex]->femaleWorms;
			double TW = (double) hostPopulation[deathIndex]->totalWorms;
			double femaleMaleWormProportion = FW/TW;

			double randNum = owner->myRandUniform();
			bool femaleWormCondition = randNum < femaleMaleWormProportion;

			// Is this worm female?
			if(femaleWormCondition)
			{
				// Remove a worm from female worms
				hostPopulation[deathIndex]->femaleWorms = hostPopulation[deathIndex]->femaleWorms - 1;
			}

			// Remove worm from total worms
			hostPopulation[deathIndex]->totalWorms = hostPopulation[deathIndex]->totalWorms - 1;
		}
	}
	else // NEW WORM
	{
		// New worm in total worms
		hostPopulation[event]->totalWorms = hostPopulation[event]->totalWorms + 1;

		double randNum = owner->myRandUniform();
		bool newFemaleWormCondition = randNum < 0.5;

		// New female worm
		if(newFemaleWormCondition)
		{
			hostPopulation[event]->femaleWorms = hostPopulation[event]->femaleWorms + 1;
		}
	}
}

// Update the freeliving populations deterministically
double CRealization::doFreeliving(double ts,double freeliving)
{
	double sumEggsOutputPerHostRho = 0.0; // Reset this value before next iteration

	for(int i=0;i<nHosts;i++)
	{
		// Female worms produce fertilised eggs only if there is a male worm around
		bool noMales = hostPopulation[i]->totalWorms == hostPopulation[i]->femaleWorms;

		productiveFemaleWorms[i] = (double) hostPopulation[i]->femaleWorms;

		if (noMales)
		{
			productiveFemaleWorms[i] = 0.0;
		}

		double lambda = owner->lambda;
		double gamma = owner->gamma;
		eggsOutputPerHost[i] = lambda*productiveFemaleWorms[i]*exp(-productiveFemaleWorms[i]*gamma);

		double rhoValue = owner->rhoValues[hostPopulation[i]->hostContactIndex];
		sumEggsOutputPerHostRho += eggsOutputPerHost[i]*rhoValue;
	}

	double psi = owner->psi;
	eggsProductionRate = (double) (psi*sumEggsOutputPerHostRho/nHosts);

	// dL/dt = K-mu*L has solution: L(0)exp(-mu*t)+K*(1-exp(-mu*t))/mu, this is exact if rate of egg production is constant in the timestep
	// Here L = freeliving, K = eggsProductionRate, mu = ReservoirDecayRate
	double ReservoirDecayRate = owner->ReservoirDecayRate;
	double expo = exp(-ReservoirDecayRate*ts);
	double freelivingNumber = freeliving*expo + (eggsProductionRate*(1.0-expo)/ReservoirDecayRate);

	return freelivingNumber;
}
