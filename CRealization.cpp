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
	surveyResultsArray = NULL;
	hostTotalWorms = NULL;
	hostFemaleWorms = NULL;
	productiveFemaleWorms = NULL;
	eggsOutputPerHost = NULL;
	eggsProductionRate = NULL;
	hostInfectionRate = NULL;
	contactAgeGroupIndex = NULL;
	treatmentAgeGroupIndex = NULL;
	rates = NULL;
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

	if (rates != NULL)
		delete[] rates;

	if(surveyResultsArray!=NULL)
	{
		for(i=0;i<owner->surveyResultTimesLength;i++)
			delete[] surveyResultsArray[i];

		delete[] surveyResultsArray;
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
	double treatmentTimeElapsed = 0;
	while((owner->treatEnd - owner->treatStart) != treatmentTimeElapsed)
	{
		localEvents.addEvent(CHEMOTHERAPY,owner->treatStart+treatmentTimeElapsed,NULL);
		treatmentTimeElapsed = treatmentTimeElapsed + owner->treatInterval;
	}
	localEvents.addEvent(CHEMOTHERAPY,owner->treatEnd,NULL); // Last treatment time

	// Add run termination point
	localEvents.addEvent(TERMINATE,owner->nYears,NULL);

	// Set up results collection for this realisation
	if(owner->surveyResultTimes!=NULL)
	{
		// There are survey results to collect
		// For each time point set up an array of surveyResultData structures the length of the population size
		surveyResultsArray = new surveyResultData*[owner->surveyResultTimesLength];
		for(int i=0;i<owner->surveyResultTimesLength;i++)
		{
			surveyResultsArray[i] = new surveyResultData[nHosts];
			// Add a pointer to the array that's going to hold the data
			localEvents.addEvent(SURVEY_EVENT,owner->surveyResultTimes[i],surveyResultsArray[i]);
		}
	}

	// DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
	localEvents.addEvent(DEBUG_EVENT,owner->nYears - 5.0,NULL);

	return true;
}

// Run the realisation
bool CRealization::run(int repNo)
{
	Event currentEvent;

	// Set some defaults
	double maxStep = 1/52; // Time steps never exceed this for deterministic update of freeliving populations

	// Some initial time-related variables
	double timeNow = 0;
	double FLlast = timeNow;
	int nextOutIndex = 0;
	double nextOutTime = owner->surveyResultTimes[nextOutIndex]; // First survey output time
	double ageingInterval = 0.25; // Checks ages every quarter year
	double nextAgeTime = ageingInterval;
	double chemoTimeElapsed = 0;
	double nextChemoTime = owner->treatStart; // First chemotherapy time
	double* nextStepCompare = new double[4];
	nextStepCompare[0] = nextOutTime;
	nextStepCompare[1] = timeNow+maxStep;
	nextStepCompare[2] = nextChemoTime;
	nextStepCompare[3] = nextAgeTime;
	double nextStep = owner->min(nextStepCompare,4); 	// Get minimum value of the nextStepCompare array


	// Delete any arrays no longer required
	if (nextStepCompare != NULL)
		delete[] nextStepCompare;

	//double elimTime = -1; // Default elimination time

	do
	{
		// Calculate sum of rates
		double sumRates = owner->sumArray(calculateRates(timeNow),nHosts+1);

		/*
		if(sumRates < 1)
		{
			elimTime = timeNow; // Next worm event is about a year away -> elimination?
		}
		 */

		double tstep = owner->myRandExponential(sumRates);

		// Do stochastic events up to the time of the next event

		if( (timeNow+tstep) < nextStep )
		{
			timeNow = timeNow + tstep;

			// Enact an event
			// Determine if dead worm or new worm
			int ratesLength = nHosts+1;
			int event = owner->multiNomBasic2(calculateRates(timeNow),ratesLength,owner->myRandUniform());

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
		else
		{
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
				int q = (int) floor(abs(timeNow - hostPopulation[i]->birthDate)/owner->demogDt);
				if (q >= owner->maxHostAge)
				{
					q = owner->maxHostAge - 1;
				}
				eggsProductionRate[i] = owner->psi*eggsOutputPerHost[i]*owner->rhoValues[owner->contactAgeGroupIndex()[q]];

				// Time step
				double ts = nextStep - FLlast;
				double expo = exp(-owner->ReservoirDecayRate*ts);
				hostPopulation[i]->freeliving = hostPopulation[i]->freeliving*expo + eggsProductionRate[i]*(1-expo)/(owner->ReservoirDecayRate);
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
					nextAgeTime = nextAgeTime + ageingInterval;
					break;

				// Apply treatment
				case CHEMOTHERAPY:
					{
						for(int i=0;i<nHosts;i++)
						{
							int q = (int) floor(abs(currentEvent.time - hostPopulation[i]->birthDate)/owner->demogDt);
							if (q >= owner->maxHostAge)
							{
								q = owner->maxHostAge - 1;
							}
							int hostContactIndex = owner->contactAgeGroupIndex()[q];
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
					}
					chemoTimeElapsed += owner->treatInterval;
					nextChemoTime = owner->treatStart + chemoTimeElapsed;
					break;

				// Any events to debug?
				case DEBUG_EVENT:
					debugEventResponse(currentEvent);
					break;

				// What to output from the simulation
				case SURVEY_EVENT:
					surveyResultResponse(currentEvent);
					nextOutIndex++;
					nextOutTime = owner->surveyResultTimes[nextOutIndex];
					break;

				// Stop simulation when maximum number of years reached
				case TERMINATE:
					break;

				default:
					owner->logStream << "Event number " << currentEvent.type << " not known.\n"
					<< "Look at #define values in CPreDetEventQueue.h for definition.\n" << std::flush;
					break;

				timeNow = nextStep;

				double* newNextStepCompare = new double[4];
				newNextStepCompare[0] = nextOutTime;
				newNextStepCompare[1] = timeNow+maxStep;
				newNextStepCompare[2] = nextChemoTime;
				newNextStepCompare[3] = nextAgeTime;
				nextStep = owner->min(newNextStepCompare,4); 	// Get new minimum value of the nextStepCompare array

				// Delete any arrays no longer required
				if (newNextStepCompare != NULL)
					delete[] newNextStepCompare;
			}
		}
	} while(currentEvent.type!=TERMINATE); // Do events until the last survey year is reached

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
	// Collect data from each host individual
	surveyResultData* outputArray = (surveyResultData*) currentEvent.subject;

	/*
	// For looking at the host ages across time
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].age = currentEvent.time - hostPopulation[i]->birthDate;
	}
	*/

	// For looking at mean female worm burdens for individual runs across time
	double sumFemaleWormsPerRun = 0;
	for(int i=0;i<nHosts;i++)
	{
		outputArray[i].nFemaleWorms += hostPopulation[i]->femaleWorms;
		sumFemaleWormsPerRun = outputArray[i].nFemaleWorms + sumFemaleWormsPerRun; // Get sum of female worms per run
	}
	outputArray->meanFemaleWormsPerRun = sumFemaleWormsPerRun/nHosts;

	return true;
}

// Just for testing...
void CRealization::debugEventResponse(Event& currentEvent)
{

}


//////////////////////////////////////////////////////////////////////////////////
/// OTHER FUNCTIONS

// Calculate the event rates
double* CRealization::calculateRates(double timeNow)
{
	// Set up some arrays
	rates = new double[nHosts+1];

	for (int i=0;i<nHosts;i++)
	{
		int q = (int) floor(abs(timeNow - hostPopulation[i]->birthDate)/owner->demogDt);
		if (q >= owner->maxHostAge)
		{
			q = owner->maxHostAge - 1;
		}
		hostInfectionRate[i] = hostPopulation[i]->freeliving*owner->myRandGamma(owner->k,owner->k)*owner->betaValues[owner->contactAgeGroupIndex()[q]];
		rates[i] = hostInfectionRate[i];
	}

	// Calculate the event rates (host infection rate and worm total death rate)
	double wormTotalDeathRate = owner->sigma*sumTotalWorms;
	rates[nHosts] = wormTotalDeathRate; // Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array

	return rates;
}
