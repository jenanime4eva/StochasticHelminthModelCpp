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
	age = 0;
	si = NULL;
	mu = NULL;
	rates = NULL;
	ratesLength = 0;
	nextStepCompare = NULL;
	nextStepCompareLength = 4;
	newNextStepCompare = NULL;
	newNextStepCompareLength = 4;
	chemoTimes = NULL;
	outTimes = NULL;
	timeNow = 0;
	surveyLength = 0;
	treatLength = 0;
	ageingInterval = 0.25; // Checks ages every quarter year
	maxStep = (double) 1/52; // Time steps never exceed this for deterministic update of freeliving populations
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

	if(si!=NULL)
		delete[] si;

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

	nHosts = owner->nHosts;

	// Set up some arrays required later
	// Total and female host worm arrays
	hostTotalWorms = new double[nHosts];
	hostFemaleWorms = new double[nHosts];
	productiveFemaleWorms = new double[nHosts];
	eggsOutputPerHost = new double[nHosts];
	hostInfectionRate = new double[nHosts];

	// Set up some variables for use later
	surveyLength = owner->surveyResultTimesLength;
	treatLength = owner->treatmentTimesLength;

	// Set up some arrays for use later
	si = new double[nHosts];
	mu = new double[nHosts];
	ratesLength = nHosts+1;
	rates = new double[ratesLength];
	nextStepCompare = new double[nextStepCompareLength];
	newNextStepCompare = new double[newNextStepCompareLength];
	chemoTimes = new double[owner->treatmentTimesLength];
	outTimes = new double[owner->surveyResultTimesLength];

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Initialise total and female worms per individual host
		// Force of infection
		si[i] = owner->myRandGamma(owner->k,owner->k); // location = 1/scale (here scale = 1/k, so location = k); shape parameter is the parameter k
		mu[i] = 2*initialWormNumber*si[i]; // These values chosen as initial conditions (initial conditions should be the equilibrium)
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
		double communityBurnIn = 1000;
		while(hostPopulation[i]->deathDate < communityBurnIn)
		{
			double newlifespan = owner->drawLifespan();
			hostPopulation[i]->birthDate = hostPopulation[i]->deathDate;
			hostPopulation[i]->deathDate = hostPopulation[i]->deathDate + newlifespan;
		}

		// Updated birth and death dates
		hostPopulation[i]->birthDate = hostPopulation[i]->birthDate - communityBurnIn;
		hostPopulation[i]->deathDate = hostPopulation[i]->deathDate - communityBurnIn;

		// Work out host ages
		hostPopulation[i]->hostAge = (int) (floor(-hostPopulation[i]->birthDate));

		// Add host death events
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);
	}

	// Add treatment events
	for(int j=0;j<treatLength;j++)
	{
		localEvents.addEvent(CHEMOTHERAPY,owner->treatmentTimes[j],NULL);
	}

	// Add run termination point
	localEvents.addEvent(TERMINATE,owner->nYears,NULL);

	// Set up results collection for this realisation
	if(owner->surveyResultTimes!=NULL)
	{
		// There are survey results to collect

		// For each time point set up an array of surveyResultData structures the length of the number of realisations
		surveyResultsArrayPerRun = new surveyResultData*[owner->surveyResultTimesLength];
		for(int i = 0; i < surveyLength; i++)
		{
			// For each time point allocate nRepetitions worth of space
			surveyResultsArrayPerRun[i] = new surveyResultData[owner->nRepetitions];
			// Add a pointer to the arrays that are going to hold the data
			localEvents.addEvent(SURVEY_EVENT,owner->surveyResultTimes[i],surveyResultsArrayPerRun[i]);
		}
	}

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,owner->nYears - 5.0,NULL);

	return true;
}

// Run the realisation
bool CRealization::run(int repNo)
{
	//printf("Run number: %d\n",repNo); // Test flag

	// Some initial time-related variables
	double timeNow = 0;
	double FLlast = timeNow;

	for(int i=0;i<surveyLength;i++)
	{
		outTimes[i] = owner->surveyResultTimes[i];
	}

	int nextOutIndex = owner->indexSmallestElement(outTimes,surveyLength);
	double nextOutTime = outTimes[nextOutIndex]; // First survey output time
	double nextAgeTime = ageingInterval;

	for(int j=0;j<treatLength;j++)
	{
		chemoTimes[j] = owner->treatmentTimes[j];
	}

	int nextChemoIndex = owner->indexSmallestElement(chemoTimes,treatLength);
	double nextChemoTime = chemoTimes[nextChemoIndex]; // First chemotherapy time

	nextStepCompare[0] = nextOutTime;
	nextStepCompare[1] = timeNow+maxStep;
	nextStepCompare[2] = nextChemoTime;
	nextStepCompare[3] = nextAgeTime;
	double nextStep = owner->min(nextStepCompare,4); // Get minimum value of the nextStepCompare array

	//double elimTime = -1; // Default elimination time

	int outCount = 0;

	Event currentEvent;
	do
	{
		// Reset these variables to zero after each iteration
		double sumTotalWorms = 0;
		double sumFemaleWorms = 0;

		// Calculate the event rates (host infection rate and worm total death rate)
		for (int i=0;i<nHosts;i++)
		{
			age = hostPopulation[i]->hostAge;
			hostInfectionRate[i] = freeliving*si[i]*owner->betaValues[owner->contactAgeGroupIndex()[age]];
			rates[i] = hostInfectionRate[i];
			sumTotalWorms += hostPopulation[i]->totalWorms; // Sum of total worms for all hosts
			sumFemaleWorms += hostPopulation[i]->femaleWorms; // Sum of female worms for all hosts
			hostTotalWorms[i] = hostPopulation[i]->totalWorms;
			hostFemaleWorms[i] = hostPopulation[i]->femaleWorms;
		}
		double wormTotalDeathRate = owner->sigma*sumTotalWorms;
		rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array

		// Calculate sum of rates
		double sumRates = owner->sumArray(rates,ratesLength);

		/*
		if(sumRates < 1)
		{
			elimTime = timeNow; // Next worm event is about a year away -> elimination?
		}
		*/

		// Calculate the time step
		double tstep = owner->myRandExponential(sumRates);

		if( (timeNow+tstep) < nextStep )
		{
			timeNow = timeNow + tstep;

			// Enact an event
			// Determine if dead worm or new worm
			int event = owner->multiNomBasic2(rates,ratesLength,owner->myRandUniform());

			if(event==ratesLength-1)
			{
				// Dead worm
				int deathIndex = owner->multiNomBasic2(hostTotalWorms,nHosts,owner->myRandUniform());

				// Are there worms present in the first place?
				if(hostPopulation[deathIndex]->totalWorms != 0)
				{
					// Is this worm female?
					if(owner->myRandUniform() < (hostFemaleWorms[deathIndex]/hostTotalWorms[deathIndex]))
					{
						hostPopulation[deathIndex]->femaleWorms = hostPopulation[deathIndex]->femaleWorms - 1; // Remove a worm
					}

					hostPopulation[deathIndex]->totalWorms = hostPopulation[deathIndex]->totalWorms - 1; // Remove a worm
				}
			}
			else
			{
				// New worm
				hostPopulation[event]->totalWorms = hostPopulation[event]->totalWorms + 1; // Add a worm

				// Female worm
				if(owner->myRandUniform() < 0.5)
				{
					hostPopulation[event]->femaleWorms = hostPopulation[event]->femaleWorms + 1; // Add a worm
				}
			}
		}
		else // (timeNow+tstep) >= nextStep
		{
			// Time step
			double ts = nextStep - FLlast;

			double sumEggsOutputPerHostRho = 0; // Reset this value before next iteration

			// Update the freeliving worm populations (worms found in the environment) deterministically
			for(int i=0;i<nHosts;i++)
			{
				productiveFemaleWorms[i] = hostPopulation[i]->femaleWorms;

				// Female worms produce fertilised eggs only if there is a male worm around
				if (hostPopulation[i]->totalWorms == hostPopulation[i]->femaleWorms)
				{
					productiveFemaleWorms[i] = 0;
				}

				eggsOutputPerHost[i] = owner->lambda*productiveFemaleWorms[i]*exp(-productiveFemaleWorms[i]*owner->gamma);

				age = hostPopulation[i]->hostAge;
				sumEggsOutputPerHostRho += eggsOutputPerHost[i]*owner->rhoValues[owner->contactAgeGroupIndex()[age]];
			}

			eggsProductionRate = owner->psi*sumEggsOutputPerHostRho/nHosts;

			// dL/dt = K-mu*L has solution: L(0)exp(-mu*t)+K*(1-exp(-mu*t))/mu, this is exact if rate of egg production is constant in the timestep
			double expo = exp(-owner->ReservoirDecayRate*ts);
			freeliving = freeliving*expo + eggsProductionRate*(1.0-expo)/(owner->ReservoirDecayRate);

			FLlast = nextStep;

			// Get next predetermined event
			localEvents.popEvent(currentEvent);

			// Do the event
			switch (currentEvent.type)
			{
				// Death of host
				case HOST_DEATH:
					hostDeathResponse(currentEvent);
					nextAgeTime = nextAgeTime + ageingInterval;
					break;

				// Apply treatment
				case CHEMOTHERAPY:
					for(int i=0;i<nHosts;i++)
					{
						// Coverage level for each host
						age = hostPopulation[i]->hostAge;
						int hostContactIndex = owner->contactAgeGroupIndex()[age];
						double individualCoverage = owner->coverage[hostContactIndex];

						// Condition for those treated, randomly chosen
						bool individualTreated = owner->myRandUniform() < individualCoverage;

						// If individualTreated condition is TRUE
						if (individualTreated)
						{
							// How many worms to die?
							double totalW = hostPopulation[i]->totalWorms;
							double femaleW = hostPopulation[i]->femaleWorms;
							double maleW = totalW - femaleW;

							double maleToDie = owner->myRandBinomial(maleW,owner->drugEff);
							double femaleToDie = owner->myRandBinomial(femaleW,owner->drugEff);

							hostPopulation[i]->totalWorms = totalW - maleToDie - femaleToDie;
							hostPopulation[i]->femaleWorms = femaleW - femaleToDie;
						}
					}
					chemoTimes[nextChemoIndex] = owner->nYears + 10; // The 10 is to make sure chemo isn't been done at maxTime
					nextChemoIndex = owner->indexSmallestElement(chemoTimes,owner->treatmentTimesLength);
					nextChemoTime = chemoTimes[nextChemoIndex];
					break;

				// Any events to debug?
				case DEBUG_EVENT:
					debugEventResponse(currentEvent);
					break;

				// What to output from the simulation
				case SURVEY_EVENT:
					surveyResultResponse(currentEvent);
					outCount++;
					outTimes[nextOutIndex] = owner->nYears + 10;
					nextOutIndex = owner->indexSmallestElement(outTimes,owner->surveyResultTimesLength);
					nextOutTime = outTimes[nextOutIndex];
					break;

				// Stop simulation when maximum number of years reached
				case TERMINATE:
					break;

				default:
					owner->logStream << "Event number " << currentEvent.type << " not known.\n"
					<< "Look at #define values in CPreDetEventQueue.h for event definition.\n" << std::flush;
					break;
			}

			timeNow = nextStep;
			newNextStepCompare[0] = nextOutTime;
			newNextStepCompare[1] = timeNow+maxStep;
			newNextStepCompare[2] = nextChemoTime;
			newNextStepCompare[3] = nextAgeTime;
			nextStep = owner->min(newNextStepCompare,4); // Get new minimum value of the nextStepCompare array
		} // End of predetermined event block
	} while((currentEvent.type!=TERMINATE) && (outCount<surveyLength)); // Do events until both nYears and the last survey year is reached

	return true;
}


//////////////////////////////////////////////////////////////////////////////////
/// PREDETERMINED EVENT RESPONSES

// Respond to host death
bool CRealization::hostDeathResponse(Event& currentEvent)
{
	// Rejuvenate host

	CHost* currentHost = (CHost*) currentEvent.subject;

	// Set birth date to now
	// Put birth slightly in the past to ensure age is just positive for catergorisation
	currentHost->birthDate = currentEvent.time - 0.001;

	// Calculate new death dates
	double lifespan = owner->drawLifespan();
	currentHost->deathDate = currentEvent.time + lifespan;

	// Calculate new force of infection (FOI)
	double newSi = owner->myRandGamma(owner->k,owner->k);
	double newMu = 2*initialWormNumber*newSi;
	currentHost->totalWorms = owner->myRandPoisson(newMu); // Total worms
	currentHost->femaleWorms = owner->myRandBinomial(currentHost->totalWorms,0.5); // Half of worm population should be female

	// Kill off all their worms
	currentHost->totalWorms = 0;
	currentHost->femaleWorms = 0;

	// New host age
	currentHost->hostAge = (int) (floor(currentEvent.time - currentHost->birthDate));

	// Assign death event
	localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

	return true;
}

// What to output from the simulation
bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// Collect data from each run
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	int runs = owner->nRepetitions; // Number of runs

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
			double currentHostAge = hostPopulation[i]->hostAge; // DEBUG: MIGHT BE CAUSING THE WORM NUMBER ISSUES

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
		outputArray[rIndex].sumFemaleWormsPerRun = sumFemaleWormsPerRun;

		// Infants
		outputArray[rIndex].sumInfantFemaleWormsPerRun = sumInfantFemaleWormsPerRun;
		outputArray[rIndex].infantNumber = infantCount;

		// Pre-SAC
		outputArray[rIndex].sumPreSACFemaleWormsPerRun = sumPreSACFemaleWormsPerRun;
		outputArray[rIndex].preSACNumber = preSACCount;

		// SAC
		outputArray[rIndex].sumSACFemaleWormsPerRun = sumSACFemaleWormsPerRun;
		outputArray[rIndex].SACNumber = SACCount;

		// Adults
		outputArray[rIndex].sumAdultFemaleWormsPerRun = sumAdultFemaleWormsPerRun;
		outputArray[rIndex].adultNumber = adultCount;
	}

	return true;
}

// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{

}

