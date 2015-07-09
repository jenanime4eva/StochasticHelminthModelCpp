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

// Class constructor
CRealization::CRealization()
{
	hostPopulation = NULL;
	owner = NULL;
	nHosts = 0;
	tinyIncrement = 0.01;
	sumTotalWorms = 0;
	sumFemaleWorms = 0;
	surveyResultsArrayPerHost = NULL;
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

	if (q != NULL)
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

	if(surveyResultsArrayPerHost!=NULL)
	{
		for(i=0;i<owner->surveyResultTimesLength;i++)
			delete[] surveyResultsArrayPerHost[i];

		delete[] surveyResultsArrayPerHost;
	}

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
	q = new int[nHosts];
	chemoTimes = new double[owner->treatmentTimesLength];
	outTimes = new double[owner->surveyResultTimesLength];

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Host lifespans
		double lifespan = owner->drawLifespan();
		hostPopulation[i]->birthDate = -owner->myRandUniform()*lifespan;
		hostPopulation[i]->deathDate = hostPopulation[i]->birthDate + lifespan;

		// Add event
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);

		// Total and female worms per individual host
		double initialWormNumber = 11.5; // Choose this number such that initial conditions are at the equilibrium
		// location = 1/scale (here scale = 1/k, so location = k); shape parameter is the parameter k
		double si = owner->myRandGamma(owner->k,owner->k);
		double mu = 2*initialWormNumber*si;
		hostPopulation[i]->totalWorms = owner->myRandPoisson(mu); // Total worms
		hostPopulation[i]->femaleWorms = owner->myRandBinomial(hostPopulation[i]->totalWorms,0.5); // Half of worm population should be female
		sumTotalWorms += hostPopulation[i]->totalWorms; // Sum of total worms for all hosts
		sumFemaleWorms += hostPopulation[i]->totalWorms; // Sum of female worms for all hosts

		// Initial freeliving worms, don't know what this number should be
		hostPopulation[i]->freeliving = 4;
	}

	// Add treatment events
	for(int i=0;i<owner->treatmentTimesLength;i++)
	{
		localEvents.addEvent(CHEMOTHERAPY,owner->treatmentTimes[i],NULL);
	}

	// Add run termination point
	localEvents.addEvent(TERMINATE,owner->nYears,NULL);

	// Set up results collection for this realisation
	if(owner->surveyResultTimes!=NULL)
	{
		// There are survey results to collect

		// For each time point set up an array of surveyResultData structures the length of the number of realisations
		surveyResultsArrayPerRun = new surveyResultData*[owner->surveyResultTimesLength];
		for(int i=0;i<owner->surveyResultTimesLength;i++)
		{
			//surveyResultsArrayPerHost[i] = new surveyResultData[nHosts];
			surveyResultsArrayPerRun[i] = new surveyResultData[owner->nRepetitions];
			// Add a pointer to the arrays that are going to hold the data
			localEvents.addEvent(SURVEY_EVENT,owner->surveyResultTimes[i],surveyResultsArrayPerRun[i]);
		}
	}

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,owner->nYears - 5.0,NULL);

	//printf ("\nrepNo = %d",repNo); // Test flag

	return true;
}

// Run the realisation
bool CRealization::run(int repNo)
{
	printf("\nrepNo: %d\n",repNo); // Test flag

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
	 // Get minimum value of the nextStepCompare array
	double nextStep = owner->min(nextStepCompare,4);

	//double elimTime = -1; // Default elimination time

	int outCount = 1;

	Event currentEvent;

	do
	{
		// Calculate the event rates (host infection rate and worm total death rate)
		for (int i=0;i<nHosts;i++)
		{
			// Work out q array
			hostAge = (int) (floor(abs(timeNow - hostPopulation[i]->birthDate))/owner->demogDt);
			q[i] = hostAge;
			if (q[i] >= owner->maxHostAge)
			{
				q[i] = owner->maxHostAge - 1;
			}

			hostInfectionRate[i] = hostPopulation[i]->freeliving*owner->myRandGamma(owner->k,owner->k)*owner->betaValues[owner->contactAgeGroupIndex()[q[i]]];
			rates[i] = hostInfectionRate[i];
		}
		double wormTotalDeathRate = owner->sigma*sumTotalWorms;
		rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array

		// Calculate sum of rates
		double sumRates = owner->sumArray(rates,nHosts+1);

		/*
		if(sumRates < 1)
		{
			elimTime = timeNow; // Next worm event is about a year away -> elimination?
		}
		 */

		double tstep = owner->myRandExponential(sumRates);

		if( (timeNow+tstep) < nextStep )
		{
			timeNow = timeNow + tstep;

			// Enact an event
			// Determine if dead worm or new worm
			int ratesLength = nHosts; // Remember in C++, array counts start from 0
			int event = owner->multiNomBasic2(rates,ratesLength,owner->myRandUniform());

			for(int i=0;i<nHosts;i++)
			{
				hostTotalWorms[i] = hostPopulation[i]->totalWorms;
				hostFemaleWorms[i] = hostPopulation[i]->femaleWorms;
			}

			if(event==ratesLength)
			{
				// Dead worm
				int deathIndex = owner->multiNomBasic2(hostTotalWorms,nHosts,owner->myRandUniform());

				// Is this worm female?
				if(owner->myRandUniform()<(hostFemaleWorms[deathIndex]/hostTotalWorms[deathIndex]))
				{
					hostPopulation[deathIndex]->femaleWorms = hostPopulation[deathIndex]->femaleWorms - 1;
				}
				hostPopulation[deathIndex]->totalWorms = hostPopulation[deathIndex]->totalWorms - 1;
			}
			else
			{
				// New worm
				hostPopulation[event]->totalWorms = hostPopulation[event]->totalWorms + 1;
				// Female worm
				if(owner->myRandUniform()<0.5)
				{
					hostPopulation[event]->femaleWorms = hostPopulation[event]->femaleWorms + 1;
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
				if ( hostPopulation[i]->totalWorms == hostPopulation[i]->femaleWorms )
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

			//double timeBarrier = nextStep + 0.001;

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
						int hostContactIndex = owner->contactAgeGroupIndex()[q[i]];
						double individualCoverage = owner->coverage[hostContactIndex];

						bool individualTreated = owner->myRandUniform()<individualCoverage;

						// If individualTreated is TRUE
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
			// Get new minimum value of the nextStepCompare array
			nextStep = owner->min(newNextStepCompare,4);
		} // End of predetermined event block
	} while((currentEvent.type!=TERMINATE) && (outCount<=surveyLength)); // Do events until both nYears and the last survey year is reached

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
	currentHost->birthDate = currentEvent.time;

	// Calculate death dates
	double lifespan = owner->drawLifespan();
	currentHost->deathDate = currentEvent.time + lifespan;

	if(currentHost->deathDate < currentEvent.time)
	{
		// Kill off all their worms
		currentHost->totalWorms = 0;
		currentHost->femaleWorms = 0;
	}

	// Assign death event
	localEvents.addEvent(HOST_DEATH,currentHost->deathDate,currentHost);

	return true;
}

// What to output from the simulation
bool CRealization::surveyResultResponse(Event& currentEvent)
{
	// Collect data from each run
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	// Set up some variables
	int reps = owner->nRepetitions;
	double sumFemaleWormsPerRun = 0;
	double sumInfantFemaleWormsPerRun = 0;
	double sumPreSACFemaleWormsPerRun = 0;
	double sumSACFemaleWormsPerRun = 0;
	double sumAdultFemaleWormsPerRun = 0;
	int infantCount = 0;
	int preSACCount = 0;
	int SACCount = 0;
	int adultCount = 0;

	// Treatment age cutoffs
	double birthCutoff = owner->treatmentBreaks[0];
	double infantCutoff = owner->treatmentBreaks[1];
	double preSACCutoff = owner->treatmentBreaks[2];
	double SACCutoff = owner->treatmentBreaks[3];
	double adultCutoff = owner->treatmentBreaks[4];

	// Female worms for each host
	for(int i=0;i<nHosts;i++)
	{
		// Looks at whole population
		sumFemaleWormsPerRun += hostPopulation[i]->femaleWorms;

		// Look at infants
		if((q[i] >= birthCutoff) && (q[i] < infantCutoff) )
		{
			sumInfantFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
			infantCount++;
		}

		// Look at pre-SAC
		if((q[i] >= infantCutoff) && (q[i] < preSACCutoff) )
		{
			sumPreSACFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
			preSACCount++;
		}

		// Look at SAC
		if((q[i] >= preSACCutoff) && (q[i] < SACCutoff) )
		{
			sumSACFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
			SACCount++;
		}

		// Look at adults
		if((q[i] >= SACCutoff) && (q[i] < adultCutoff) )
		{
			sumAdultFemaleWormsPerRun += hostPopulation[i]->femaleWorms;
			adultCount++;
		}

	}

	// For each realisation
	for(int rIndex=0;rIndex<reps;rIndex++)
	{
		outputArray[rIndex].meanFemaleWormsPerRun = (double) sumFemaleWormsPerRun/nHosts;

		if(infantCount!=0)
		{
			outputArray[rIndex].meanInfantFemaleWormsPerRun = (double) sumInfantFemaleWormsPerRun/infantCount;
		}
		else
		{
			outputArray[rIndex].meanInfantFemaleWormsPerRun = 0;
		}

		if (preSACCount!=0)
		{
			outputArray[rIndex].meanPreSACFemaleWormsPerRun = (double) sumPreSACFemaleWormsPerRun/preSACCount;
		}
		else
		{
			outputArray[rIndex].meanPreSACFemaleWormsPerRun = 0;
		}

		if (SACCount!=0)
		{
			outputArray[rIndex].meanSACFemaleWormsPerRun = (double) sumSACFemaleWormsPerRun/SACCount;
		}
		else
		{
			outputArray[rIndex].meanSACFemaleWormsPerRun = 0;
		}

		if (adultCount!=0)
		{
			outputArray[rIndex].meanAdultFemaleWormsPerRun = (double) sumAdultFemaleWormsPerRun/adultCount;
		}
		else
		{
			outputArray[rIndex].meanAdultFemaleWormsPerRun = 0;
		}
	}

	return true;
}

// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{

}

