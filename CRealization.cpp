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
	freeliving = 4.0; // Initial freeliving worms, don't know what this number should be
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
	timeNow = 0;
	birthCutoff = 0;
	infantCutoff = 0;
	preSACCutoff = 0;
	SACCutoff = 0;
	adultCutoff = 0;
	counter = 0; // Global variable counter for testing purposes
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
	productiveFemaleWorms = new int[nHosts];
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
		hostPopulation[i]->si = owner->myRandGamma(owner->k,owner->k); // location = 1/scale (here scale = 1/k, so location = k); shape parameter is the parameter k
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
		hostPopulation[i]->hostAge = -hostPopulation[i]->birthDate; // 0 - birthDate (as time at this point is zero)

		int age = (int) floor(hostPopulation[i]->hostAge);

		// Work out hostContactIndex and hostTreatIndex
		hostPopulation[i]->hostContactIndex = owner->contactAgeGroupIndex()[age];
		hostPopulation[i]->hostTreatIndex = owner->treatmentAgeGroupIndex()[age];
	}

	// Set up results collection for this realization.
	if(owner->surveyResultTimes!=NULL)
	{
		// there are survey results to collect.
		// for each time point set up an array of surveyResultData structures the length of the population size
		surveyResultsArrayPerRun = new surveyResultData*[owner->surveyResultTimesLength];
		for(int i=0;i<owner->surveyResultTimesLength;i++)
		{
			surveyResultsArrayPerRun[i] = new surveyResultData[nHosts];
			localEvents.addEvent(SURVEY_EVENT,owner->surveyResultTimes[i],surveyResultsArrayPerRun[i]); // Add a pointer to the array that's going to hold the data
		}
	}

	return true;
}

// Run the realisation
bool CRealization::run(int repNo)
{
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

	// Get minimum value of compareArray
	compareArray[0] = nextOutTime;
	compareArray[1] = timeNow+maxStep;
	compareArray[2] = nextChemoTime;
	compareArray[3] = nextAgeTime;
	double nextStep = owner->min(compareArray,4);

	int outCount = 0;

	do
	{
		// Store individual host total worms into the hostTotalWorms array
		for (int i=0;i<nHosts;i++)
		{
			hostTotalWorms[i] = (double) hostPopulation[i]->totalWorms;
		}

		// Work out sum of the total worms array
		int sumTotalWorms = (int) owner->sumArray(hostTotalWorms,nHosts);

		// Calculate the event rates (host infection rate and worm total death rate)
		calculateEventRates(sumTotalWorms);

		// Calculate sum of rates
		double sumRates = owner->sumArray(rates,ratesLength);

		// Calculate the time step
		double tstep = owner->myRandExponential(sumRates);
		//printf("tstep: %f\n",tstep);

		if( (timeNow+tstep) < nextStep )
		{
			timeNow = timeNow + tstep;

			// Enact an event
			doEvent(hostTotalWorms);
		}
		else // (timeNow+tstep) >= nextStep
		{
			// Time step
			ts = nextStep - FLlast;
			//printf("ts: %f\n",ts);

			// Update the freeliving worm populations (worms found in the environment) deterministically
			freeliving = doFreeliving(ts,freeliving);

			FLlast = nextStep;

			double timeBarrier = nextStep + 0.001;

			// Ageing and death
			if(timeBarrier > nextAgeTime)
			{
				doDeath(timeNow);
				nextAgeTime = nextAgeTime + ageingInterval;
			}

			// Chemotherapy
			if(timeBarrier > nextChemoTime)
			{
				doChemo();
				chemoTimes[nextChemoIndex] = nYears + 10; // The 10 is to make sure chemo isn't been done at maxTime
				nextChemoIndex = owner->indexSmallestElement(chemoTimes,treatLength);
				nextChemoTime = chemoTimes[nextChemoIndex];
			}

			// Record
			if(timeBarrier > nextOutTime)
			{
				outCount++; // Add one to outCount
				outTimes[nextOutIndex] = nYears + 10;
				nextOutIndex = owner->indexSmallestElement(outTimes,surveyLength);
				nextOutTime = outTimes[nextOutIndex];
			}

			// Get next predetermined event
			localEvents.popEvent(currentEvent);

			// Do the event
			switch (currentEvent.type)
			{
				// What to output from the simulation
				case SURVEY_EVENT:
					//printf("survey currentEvent.time: %f\n",currentEvent.time);
					surveyResultResponse(currentEvent);
					break;

				/*
				default:
					owner->logStream << "Event number " << currentEvent.type << " not known.\n"
					<< "Look at #define values in CPreDetEventQueue.h for event definition.\n" << std::flush;
					break;
				*/
			}

			timeNow = nextStep;

			// Get new minimum of compareArray
			compareArray[0] = nextOutTime;
			compareArray[1] = timeNow+maxStep;
			compareArray[2] = nextChemoTime;
			compareArray[3] = nextAgeTime;
			nextStep = owner->min(compareArray,4);

			// TODO
			//printf("timeNow: %f\n",timeNow);
			//printf("currentEvent.time: %f\n",currentEvent.time);
			//printf("nextStep %f\n",nextStep);
			//printf("hostPopulation[0] %d\n",hostPopulation[0]->femaleWorms);

		} // End of predetermined event block
	} while((timeNow<nYears) && (outCount<surveyLength)); // Do events until both nYears and the last survey year is reached

	return true;
}


//////////////////////////////////////////////////////////////////////////////////
/// PREDETERMINED EVENT RESPONSES

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


//////////////////////////////////////////////////////////////////////////////////
/// OTHER FUNCTIONS

// Calculate the event rates
void CRealization::calculateEventRates(int sumTotalWorms)
{
	for (int i=0;i<nHosts;i++)
	{
		hostInfectionRate[i] = freeliving*hostPopulation[i]->si*owner->betaValues[hostPopulation[i]->hostContactIndex];
		rates[i] = hostInfectionRate[i];
	}
	double wormTotalDeathRate = (double) owner->sigma*sumTotalWorms;
	rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array
}

// Enact an event
void CRealization::doEvent(double* hostTotalWorms)
{
	// Determine if dead worm or new worm
	int event = owner->multiNomBasic2(rates,ratesLength);
	//printf("event %d\n",event);
	if(event==ratesLength-1) // The -1 is there because C++ counts array entries from 0
	{
		// Dead worm
		int deathIndex = owner->multiNomBasic2(hostTotalWorms,nHosts);

		double femaleMaleWormProportion = (double) (hostPopulation[deathIndex]->femaleWorms / hostPopulation[deathIndex]->totalWorms);

		bool femaleWormCondition = owner->myRandUniform() < femaleMaleWormProportion;

		// Is this worm female?
		if(femaleWormCondition)
		{
			// Remove a worm from female worms
			hostPopulation[deathIndex]->femaleWorms = hostPopulation[deathIndex]->femaleWorms - 1;
		}

		// Remove worm from total worms
		hostPopulation[deathIndex]->totalWorms = hostPopulation[deathIndex]->totalWorms - 1;
	}
	else
	{
		// New worm in total worms
		hostPopulation[event]->totalWorms = hostPopulation[event]->totalWorms + 1;

		bool newFemaleWormCondition = owner->myRandUniform() < 0.5;

		// New female worm
		if(newFemaleWormCondition)
		{
			hostPopulation[event]->femaleWorms = hostPopulation[event]->femaleWorms + 1;
		}
	}
}

// Update the freeliving populations deterministically
double CRealization::doFreeliving(double ts,int freeliving)
{
	double sumEggsOutputPerHostRho = 0.0; // Reset this value before next iteration

	for(int i=0;i<nHosts;i++)
	{
		// Female worms produce fertilised eggs only if there is a male worm around
		bool noMales = hostPopulation[i]->totalWorms == hostPopulation[i]->femaleWorms;

		productiveFemaleWorms[i] = hostPopulation[i]->femaleWorms;

		if (noMales)
		{
			productiveFemaleWorms[i] = 0;
		}

		eggsOutputPerHost[i] = (double) owner->lambda*productiveFemaleWorms[i]*exp(-productiveFemaleWorms[i]*owner->gamma);

		sumEggsOutputPerHostRho += eggsOutputPerHost[i]*owner->rhoValues[hostPopulation[i]->hostContactIndex];
	}

	eggsProductionRate = owner->psi*sumEggsOutputPerHostRho/nHosts;

	// dL/dt = K-mu*L has solution: L(0)exp(-mu*t)+K*(1-exp(-mu*t))/mu, this is exact if rate of egg production is constant in the timestep
	// Here L = freeliving, K = eggsProductionRate, mu = ReservoirDecayRate
	double expo = exp(-owner->ReservoirDecayRate*ts);
	double freelivingNumber = freeliving*expo + eggsProductionRate*(1.0-expo)/(owner->ReservoirDecayRate);

	return freelivingNumber;
}

// Respond to host death
void CRealization::doDeath(double timeNow)
{
	for(int i=0;i<nHosts;i++)
	{
		// Identify the indices of the dead
		bool theDead = hostPopulation[i]->deathDate < timeNow;

		if(theDead)
		{
			// Set birth date to now
			// Put birth slightly in the past to ensure age is just positive for catergorisation
			hostPopulation[i]->birthDate = timeNow - 0.001;

			// Calculate new death dates
			double lifespan = owner->drawLifespan();
			hostPopulation[i]->deathDate = timeNow + lifespan;

			// Calculate new force of infection (FOI)
			hostPopulation[i]->si = owner->myRandGamma(owner->k,owner->k);

			// Kill off all their worms
			hostPopulation[i]->totalWorms = 0;
			hostPopulation[i]->femaleWorms = 0;

			// New host age
			hostPopulation[i]->hostAge = timeNow - hostPopulation[i]->birthDate;

			int age = (int) floor(hostPopulation[i]->hostAge);

			// Update hostContactIndex and hostTreatIndex
			hostPopulation[i]->hostContactIndex = owner->contactAgeGroupIndex()[age];
			hostPopulation[i]->hostTreatIndex = owner->treatmentAgeGroupIndex()[age];
		}
	}
}

// Chemotherapy
void CRealization::doChemo()
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
			int totalW = hostPopulation[i]->totalWorms;
			int femaleW = hostPopulation[i]->femaleWorms;
			int maleW = totalW - femaleW;

			double drugEfficacy = owner->drugEff;

			int maleToDie = owner->myRandBinomial(maleW,drugEfficacy);
			int femaleToDie = owner->myRandBinomial(femaleW,drugEfficacy);

			hostPopulation[i]->totalWorms = hostPopulation[i]->totalWorms - maleToDie - femaleToDie;
			hostPopulation[i]->femaleWorms = hostPopulation[i]->femaleWorms - femaleToDie;
		}
	}
}
