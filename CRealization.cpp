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
	sumTotalWorms = 0;
	sumFemaleWorms = 0;
	surveyResultsArrayPerRun = NULL;
	hostTotalWorms = NULL;
	hostFemaleWorms = NULL;
	productiveFemaleWorms = NULL;
	eggsOutputPerHost = NULL;
	eggsProductionRate = NULL;
	hostInfectionRate = NULL;
	contactAgeGroupIndex = NULL;
	treatmentAgeGroupIndex = NULL;
	q = NULL;
	hostAge = 0;
	rates = NULL;
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

	if (eggsProductionRate != NULL)
		delete[] eggsProductionRate;

	if (hostInfectionRate != NULL)
		delete[] hostInfectionRate;

	if (contactAgeGroupIndex != NULL)
		delete[] contactAgeGroupIndex;

	if (treatmentAgeGroupIndex != NULL)
		delete[] treatmentAgeGroupIndex;

	if(q!=NULL)
		delete[] q;


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
	eggsProductionRate = new double[nHosts];
	hostInfectionRate = new double[nHosts];

	// Set up some variables for use later
	surveyLength = owner->surveyResultTimesLength;
	treatLength = owner->treatmentTimesLength;

	// Set up some arrays for use later
	rates = new double[nHosts+1];
	nextStepCompare = new double[nextStepCompareLength];
	newNextStepCompare = new double[newNextStepCompareLength];
	q = new int[nHosts]; // To store array of host ages
	chemoTimes = new double[owner->treatmentTimesLength];
	outTimes = new double[owner->surveyResultTimesLength];

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Initialise total and female worms per individual host
		double si = owner->myRandGamma(owner->k,owner->k); // location = 1/scale (here scale = 1/k, so location = k); shape parameter is the parameter k
		double mu = 2*initialWormNumber*si; // These values chosen as initial conditions (initial conditions should be the equilibrium)
		hostPopulation[i]->totalWorms = owner->myRandPoisson(mu); // Total worms
		hostPopulation[i]->femaleWorms = owner->myRandBinomial(hostPopulation[i]->totalWorms,0.5); // Half of worm population should be female
		sumTotalWorms += hostPopulation[i]->totalWorms; // Sum of total worms for all hosts
		sumFemaleWorms += hostPopulation[i]->femaleWorms; // Sum of female worms for all hosts

		// Initial freeliving worms, don't know what this number should be
		hostPopulation[i]->freeliving = 4.0;

		// Distribute birth dates such that at time zero have the right distribution of ages
		// Give a sample of ages at the start of the simulation
		// For an exponential distribution, the survival function S(t) = exp(-t/mu)
		// So to sample from this distribution, generate a random number on 0 1 and then invert this function

		// Host life span
		double lifespan = owner->drawLifespan();

		// Host birth date
		hostPopulation[i]->birthDate = -owner->myRandUniform()*lifespan;

		// Host death date
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;

		// Equilibrate the population first, in the absence of understanding how to generate it in the first place
		double communityBurnIn = 1000;
		while(hostPopulation[i]->deathDate < communityBurnIn)
		{
			hostPopulation[i]->birthDate = hostPopulation[i]->deathDate;
			double newlifespan = owner->drawLifespan();
			hostPopulation[i]->deathDate = hostPopulation[i]->deathDate + newlifespan;
		}

		// Updated birth and death dates
		hostPopulation[i]->birthDate = hostPopulation[i]->birthDate - communityBurnIn;
		hostPopulation[i]->deathDate = hostPopulation[i]->deathDate - communityBurnIn;

		// Work out q array of host ages
		hostAge = (int) (floor(-hostPopulation[i]->birthDate));
		q[i] = hostAge;

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
		// Calculate the event rates (host infection rate and worm total death rate)
		for (int i=0;i<nHosts;i++)
		{
			hostInfectionRate[i] = hostPopulation[i]->freeliving*owner->myRandGamma(owner->k,owner->k)*owner->betaValues[owner->contactAgeGroupIndex()[q[i]]];
			rates[i] = hostInfectionRate[i];
		}
		double wormTotalDeathRate = owner->sigma*sumTotalWorms;
		rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array

		// Calculate sum of rates
		int ratesLength = nHosts; // Remember in C++, array counts start from 0 so no need for a "+1" here
		double sumRates = owner->sumArray(rates,ratesLength);

		/*
		if(sumRates < 1)
		{
			elimTime = timeNow; // Next worm event is about a year away -> elimination?
		}
		*/

		for (int i=0;i<nHosts;i++)
		{
			hostTotalWorms[i] = hostPopulation[i]->totalWorms;
			hostFemaleWorms[i] = hostPopulation[i]->femaleWorms;
			//printf("%f\n",hostPopulation[i]->femaleWorms);
		}

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
					if(owner->myRandUniform()<(hostFemaleWorms[deathIndex]/hostTotalWorms[deathIndex]))
					{
						hostPopulation[deathIndex]->femaleWorms = hostPopulation[deathIndex]->femaleWorms - 1; // Remove a worm
					}
					hostPopulation[deathIndex]->totalWorms = hostPopulation[deathIndex]->totalWorms - 1; // Remove a worm
				}
				else // There were no worms left to begin with
				{
					hostPopulation[deathIndex]->totalWorms = 0;
					hostPopulation[deathIndex]->femaleWorms = 0;
				}
			}
			else
			{
				// New worm
				hostPopulation[event]->totalWorms = hostPopulation[event]->totalWorms + 1; // Add a worm
				// Female worm
				if(owner->myRandUniform()<0.5)
				{
					hostPopulation[event]->femaleWorms = hostPopulation[event]->femaleWorms + 1; // Add a worm
				}
			}
		}
		else // (timeNow+tstep) >= nextStep
		{
			// Time step
			double ts = nextStep - FLlast;

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
				eggsProductionRate[i] = owner->psi*eggsOutputPerHost[i]*owner->rhoValues[owner->contactAgeGroupIndex()[q[i]]];
				double expo = exp(-owner->ReservoirDecayRate*ts);
				hostPopulation[i]->freeliving = hostPopulation[i]->freeliving*expo + eggsProductionRate[i]*(1.0-expo)/(owner->ReservoirDecayRate);
			}

			FLlast = nextStep;

			// Get next predetermined event
			localEvents.popEvent(currentEvent);

			// Do the event
			switch (currentEvent.type)
			{
				// Death of host
				case HOST_DEATH:
					hostDeathResponse(currentEvent);
					for(int i=0;i<nHosts;i++)
					{
						// Work out new q array of host ages (to update contact age categories)
						hostAge = (int) (floor(currentEvent.time - hostPopulation[i]->birthDate));
						q[i] = hostAge;
					}
					nextAgeTime = nextAgeTime + ageingInterval;
					break;

				// Apply treatment
				case CHEMOTHERAPY:
					for(int i=0;i<nHosts;i++)
					{
						int hostContactIndex = owner->contactAgeGroupIndex()[q[i]];
						double individualCoverage = owner->coverage[hostContactIndex];

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
					<< "Look at #define values in CPreDetEventQueue.h for definition.\n" << std::flush;
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

	// Calculate new force of infections (FOIs)
	double si = owner->myRandGamma(owner->k,owner->k);
	double mu = 2*initialWormNumber*si;
	currentHost->totalWorms = owner->myRandPoisson(mu); // Total worms
	currentHost->femaleWorms = owner->myRandBinomial(currentHost->totalWorms,0.5); // Half of worm population should be female

	// Kill off all their worms
	currentHost->totalWorms = 0;
	currentHost->femaleWorms = 0;

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

		// Set up current age of host variable
		double currentHostAge = 0;

		// Female worms for each host
		for(int i=0;i<nHosts;i++)
		{
			// Get current age of host
			currentHostAge = currentEvent.time - hostPopulation[i]->birthDate;

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
		//printf("sumPreSACFemaleWormsPerRun %f\n",sumPreSACFemaleWormsPerRun);

		// Whole population
		outputArray[rIndex].meanFemaleWormsPerRun = (double) sumFemaleWormsPerRun/nHosts;

		// Infants
		outputArray[rIndex].meanInfantFemaleWormsPerRun = (double) sumInfantFemaleWormsPerRun/(infantCount + 0.01); // 0.01 is to avoid division by zero

		// Pre-SAC
		outputArray[rIndex].meanPreSACFemaleWormsPerRun = (double) sumPreSACFemaleWormsPerRun/(preSACCount + 0.01); // 0.01 is to avoid division by zero

		// SAC
		outputArray[rIndex].meanSACFemaleWormsPerRun = (double) sumSACFemaleWormsPerRun/(SACCount + 0.01); // 0.01 is to avoid division by zero

		// Adults
		outputArray[rIndex].meanAdultFemaleWormsPerRun = (double) sumAdultFemaleWormsPerRun/(adultCount + 0.01); // 0.01 is to avoid division by zero
	}

	return true;
}

// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{

}

