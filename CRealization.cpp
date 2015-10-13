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
	lifespan = 0;
	tinyIncrement = 0.01;
	freeliving = 0;
	compareArray = NULL;
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
	chemoTimes = NULL;
	outTimes = NULL;
	surveyLength = 0;
	treatLength = 0;
	birthCutoff = 0;
	infantCutoff = 0;
	preSACCutoff = 0;
	SACCutoff = 0;
	adultCutoff = 0;
	treatmentOnOff = 0;
	treatStart = 0;
	treatEnd = 0;
	treatInterval = 0;
	sigma = 0;
	psi = 0;
	k = 0;
	drugEfficacy = 0;
	lambda = 0;
	gamma = 0;
	ReservoirDecayRate = 0;
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

	if (compareArray != NULL)
		delete[] compareArray;

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
	treatmentOnOff = owner->treatmentOnOff; // Is treatment on (value is 1) or off (value is 0)?
	treatStart = owner->treatStart; // Treatment start year
	treatEnd = owner->treatEnd; // Treatment end year
	treatInterval= owner->treatInterval; // Treatment interval
	sigma = owner->sigma; // Worm death rate
	psi = owner->psi;
	k = owner->k;
	drugEfficacy = owner->drugEff; // Drug efficacy
	lambda = owner->lambda; // EPG
	gamma = owner->gamma; // Exponential density dependence of parasite adult stage
	ReservoirDecayRate = owner->ReservoirDecayRate; // decay rate of eggs in the environment

	// Set up some arrays required later
	compareArray = new double[4];
	hostTotalWorms = new double[nHosts];
	productiveFemaleWorms = new int[nHosts];
	eggsOutputPerHost = new double[nHosts];
	hostInfectionRate = new double[nHosts];

	// Set up some arrays for use later
	mu = new double[nHosts];
	ratesLength = nHosts+1;
	rates = new double[ratesLength];
	chemoTimes = new double[treatLength];
	outTimes = new double[surveyLength];

	// Treatment age cutoffs
	birthCutoff = owner->treatmentBreaks[0];
	infantCutoff = owner->treatmentBreaks[1];
	preSACCutoff = owner->treatmentBreaks[2];
	SACCutoff = owner->treatmentBreaks[3];
	adultCutoff = owner->treatmentBreaks[4];

	double initialWormNumber = 11.5; // Choose this number such that initial conditions are at the equilibrium

	// Set up host population array
	hostPopulation = new CHost* [nHosts];
	for (int i=0;i<nHosts;i++)
	{
		hostPopulation[i] = new CHost;

		// Initialise total and female worms per individual host
		// Force of infection
		double siValue = owner->myRandGamma(k,k); // location = 1/scale (here scale = 1/k, so location = k); shape parameter is the parameter k
		hostPopulation[i]->si = siValue;
		mu[i] = 2.0*initialWormNumber*hostPopulation[i]->si; // These values chosen as initial conditions (initial conditions should be the equilibrium)
		int TW = owner->myRandPoisson(mu[i]);
		hostPopulation[i]->totalWorms = TW; // Total worms
		int FW = owner->myRandBinomial(TW,0.5);
		hostPopulation[i]->femaleWorms = FW; // Half of worm population should be female

		// Distribute birth dates such that at time zero have the right distribution of ages
		// Give a sample of ages at the start of the simulation
		// For an exponential distribution, the survival function S(t) = exp(-t/mu)
		// So to sample from this distribution, generate a random number on 0 1 and then invert this function

		// Host life span
		lifespan = owner->drawLifespan();

		// Host birth and death date
		double randNum = owner->myRandUniform();
		double birthDate = -randNum*lifespan;
		double deathDate = birthDate + lifespan;

		// Equilibrate the population first, in the absence of understanding how to generate it in the first place
		double communityBurnIn = 1000.0; // 1000 years should be plenty of time to equilibrate the population

		double currentBirthDate = birthDate;
		double currentDeathDate = deathDate;

		while(currentDeathDate < communityBurnIn)
		{
			lifespan = owner->drawLifespan();
			currentBirthDate = currentDeathDate;
			currentDeathDate = currentBirthDate + lifespan;
		}

		// Updated birth and death dates
		hostPopulation[i]->birthDate = currentBirthDate - communityBurnIn;
		hostPopulation[i]->deathDate = currentDeathDate - communityBurnIn;

		// Host age
		double hostAge = -hostPopulation[i]->birthDate; // 0 - birthDate (as time at this point is zero)

		// Work out hostContactIndex and hostTreatIndex
		int age = (int) ceil(hostAge);
		int contactIndex =  owner->contactAgeGroupIndex()[age];
		int treatIndex =  owner->treatmentAgeGroupIndex()[age];
		hostPopulation[i]->hostContactIndex = contactIndex;
		hostPopulation[i]->hostTreatIndex = treatIndex;

		// Add host death events
		localEvents.addEvent(HOST_DEATH,hostPopulation[i]->deathDate,hostPopulation[i]);
	}

	freeliving = initialWormNumber; // Initial freeliving worms, don't know what this number should be

	for(int j=0;j<treatLength;j++)
	{
		double currentTreatmentTime = owner->treatmentTimes[j];
		chemoTimes[j] = currentTreatmentTime;
	}

	// Add treatment events
	for(int j=0;j<treatLength;j++)
	{
		// Add treatment events just slightly before the specified treat times to make sure treatment has definitely been carried out by those times
		localEvents.addEvent(CHEMOTHERAPY,chemoTimes[j] - 0.001,NULL);
	}

	for(int i=0;i<surveyLength;i++)
	{
		double currentOutTime = owner->surveyResultTimes[i];
		outTimes[i] = currentOutTime;
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
	// Print run number
	//printf("RUN NUMBER: %d\n",repNo);

	// Time right now
	double timeNow = 0;

	// Time between reservoir update
	double deltaTimeDet = (double) 1/52; // Time steps never exceed this for deterministic update of freeliving populations

	 // Time of the next reservoir update
	double timeRes = timeNow + deltaTimeDet;

	bool reservoirNext; // This variable will be either TRUE (value of 1) or FALSE (value of 0)

	Event nextEvent; // Take nextEvent off queue and store it here

	// Allocate values to compareArray
	compareArray[0] = outTimes[0]; // First survey time
	compareArray[1] = 0.25; // Ageing interval
	compareArray[2] = chemoTimes[0]; // First treat time
	compareArray[3] = timeRes; // Time of next reservoir update

	nextEvent.time = owner->min(compareArray,4); // For each run, reset this value to the minimum of compareArray

	// Which event?
	double timeEvent = timeRes;
	reservoirNext = 1; // reservoirNext is TRUE
	if(timeRes > nextEvent.time)
	{
		timeEvent = nextEvent.time;
		reservoirNext = 0; // reservoirNext is FALSE
	}

	do
	{
		// Update hostContactIndex and hostTreatIndex
		for(int i=0;i<nHosts;i++)
		{
			// Host age
			double hostAge = timeNow - hostPopulation[i]->birthDate;

			// Work out hostContactIndex and hostTreatIndex
			int age = (int) ceil(hostAge);
			int contactIndex =  owner->contactAgeGroupIndex()[age];
			int treatIndex =  owner->treatmentAgeGroupIndex()[age];
			hostPopulation[i]->hostContactIndex = contactIndex;
			hostPopulation[i]->hostTreatIndex = treatIndex;
		}

		// Calculate the event rates (host infection rate and worm total death rate)
		calculateEventRates();

		// Last rates array entry
		double lastRate = rates[ratesLength-1];
		//printf("timeNow: %f	hostInfectionRate: %f	wormTotalDeathRate: %f\n",timeNow,rates[nHosts-1],rates[nHosts]-rates[nHosts-1]);

		// Stochastic time step
		double deltaTimeStoch = owner->myRandExponential(lastRate);

		// Time of next event
		double timeNextEvent = timeNow + deltaTimeStoch;

		while(timeNextEvent < timeEvent)
		{
			// Enact a worm event
			doEvent();

			// Calculate rates and update deltaTimeStoch
			calculateEventRates();
			lastRate = rates[ratesLength-1];
			deltaTimeStoch = owner->myRandExponential(lastRate);
			timeNow = timeNextEvent;
			timeNextEvent = timeNow + deltaTimeStoch;
		}

		// timeNextEvent must be greater or equal to timeEvent now
		// Time for the event

		// If reservoirNext is TRUE
		if(reservoirNext)
		{
			// It's deltaTimeDet since the last reservoir update, so it's time to update the reservoir

			// Update the freeliving worm populations (worms found in the environment) deterministically
			freeliving = doFreeliving(timeNow,freeliving); // TODO: CHECK THIS
			//printf("timeNow: %f freeliving: %f\n",timeNow,freeliving);

			timeNow = timeRes;
			timeRes = timeNow + deltaTimeDet;
		}
		else // reservoirNext is FALSE
		{
			// Time for the queued events to take place
			timeNow = nextEvent.time;

			// Do the event
			switch (nextEvent.type)
			{
				// Death of host
				case HOST_DEATH:
					hostDeathResponse(nextEvent);
					break;

				// Apply treatment
				case CHEMOTHERAPY:
					//printf("BEFORE %d\t",hostPopulation[10]->femaleWorms);
					hostChemoResponse(nextEvent);
					//printf("AFTER %d\n",hostPopulation[10]->femaleWorms);
					break;

				// Any events to debug?
				case DEBUG_EVENT:
					debugEventResponse(nextEvent);
					break;

				// What to output from the simulation
				case SURVEY_EVENT:
					surveyResultResponse(nextEvent);
					break;

				default:
					owner->logStream << "Event number " << nextEvent.type << " not known.\n"
					<< "Look at #define values in CPreDetEventQueue.h for event definition.\n" << std::flush;
					break;
			}

			// Update nextEvent
			// Get next predetermined event
			localEvents.popEvent(nextEvent);

		} // End of predetermined event block

		// Which event?
		timeEvent = timeRes;
		reservoirNext = 1; // reservoirNext is TRUE
		if(timeRes > nextEvent.time)
		{
			timeEvent = nextEvent.time;
			reservoirNext = 0; // reservoirNext is FALSE
		}

		//printf("timeNow: %f\n",timeNow);
		//printf("nextEvent.time: %f\n",nextEvent.time);

	} while(timeNow <= nYears); // Do events until nYears is reached

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
	lifespan = owner->drawLifespan();

	// Set birth date to now
	// Put birth slightly in the past to ensure age is just positive for categorisation
	currentHost->birthDate = currentEvent.time - 0.001;

	// Calculate new death dates
	currentHost->deathDate = currentEvent.time + lifespan;

	// Calculate new force of infection (FOI)
	double siValue = owner->myRandGamma(k,k);
	currentHost->si = siValue;

	// Kill off all their worms
	currentHost->totalWorms = 0;
	currentHost->femaleWorms = 0;

	// Update hostContactIndex and hostTreatIndex
	double hostAge = currentEvent.time - currentHost->birthDate;
	int age = (int) ceil(hostAge);
	int contactIndex = owner->contactAgeGroupIndex()[age];
	int treatIndex = owner->treatmentAgeGroupIndex()[age];
	currentHost->hostContactIndex = contactIndex;
	currentHost->hostTreatIndex = treatIndex;

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
		int treatIndex = hostPopulation[i]->hostTreatIndex;
		double individualCoverage = owner->coverage[treatIndex];

		// Condition for those treated, randomly chosen
		double randNum = owner->myRandUniform();
		bool individualTreated = randNum < individualCoverage;

		// If individualTreated condition is TRUE
		if (individualTreated)
		{
			// How many worms to die?
			int TW = hostPopulation[i]->totalWorms;
			int FW = hostPopulation[i]->femaleWorms;
			int MW = TW - FW;

			int maleToDie = owner->myRandBinomial(MW,drugEfficacy);
			int femaleToDie = owner->myRandBinomial(FW,drugEfficacy);

			hostPopulation[i]->totalWorms = TW - maleToDie - femaleToDie;
			hostPopulation[i]->femaleWorms = FW - femaleToDie;
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
		double currentHostAge = outputArray[i].hostAge;

		// Look at whole population
		outputArray[i].femaleWorms = hostPopulation[i]->femaleWorms; // Number of female worms for ith host
		outputArray[i].totalHostNumber = 1;

		// Look at infants
		bool infants = (currentHostAge > birthCutoff) && (currentHostAge <= infantCutoff);

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
		bool preSAC = (currentHostAge > infantCutoff) && (currentHostAge <= preSACCutoff);

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
		bool SAC = (currentHostAge > preSACCutoff) && (currentHostAge <= SACCutoff);

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
		bool adults = (currentHostAge > SACCutoff) && (currentHostAge <= adultCutoff);

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
		int contactIndex =  hostPopulation[i]->hostContactIndex;
		double betaValue = owner->betaValues[contactIndex];
		hostInfectionRate[i] = freeliving*hostPopulation[i]->si*betaValue;
		rates[i] = cumul + hostInfectionRate[i];
		cumul = rates[i];
	}

	double wormTotalDeathRate = (double) (sigma*sumTotalWorms);

	// Append wormTotalDeathRate onto the hostInfectionRate array to complete rates array
	rates[nHosts] = cumul + wormTotalDeathRate;
}

// Enact an event
void CRealization::doEvent()
{
	// Determine if dead worm or new worm
	// If event is equal to 0,1,2,...,nHosts-1 it's a NEW WORM, otherwise (if event is equal to nHost) it's a DEAD WORM
	double currentRand = owner->myRandUniform();
	int event = owner->multiNomBasic(rates,ratesLength,currentRand);
	//printf("event=%d\n",event);

	if(event==nHosts) // DEAD WORM
	{
		// Now let's make hostTotalWorms a cumulative array (so that multiNomBasic may be used on it in a bit)
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
double CRealization::doFreeliving(double timeNow,double freeliving)
{
	double sumEggsOutputPerHostRho = 0.0; // Reset this value before next iteration

	for(int i=0;i<nHosts;i++)
	{
		// Female worms produce fertilised eggs only if there is a male worm
		bool noMales = hostPopulation[i]->totalWorms == hostPopulation[i]->femaleWorms;

		productiveFemaleWorms[i] = hostPopulation[i]->femaleWorms;

		if (noMales)
		{
			productiveFemaleWorms[i] = 0;
		}

		eggsOutputPerHost[i] = (double) lambda*productiveFemaleWorms[i]*exp(-productiveFemaleWorms[i]*gamma);

		int contactIndex = hostPopulation[i]->hostContactIndex;
		double rhoValue = owner->rhoValues[contactIndex];
		sumEggsOutputPerHostRho += eggsOutputPerHost[i]*rhoValue;
	}

	eggsProductionRate = (double) (psi*sumEggsOutputPerHostRho/nHosts);

	// dL/dt = K-mu*L has solution L=L(0)exp(-mu*t)+K*(1-exp(-mu*t))/mu, this is exact if rate of egg production is constant in the timestep
	// Here L = freeliving, K = eggsProductionRate, mu = ReservoirDecayRate
	double expo = exp(-ReservoirDecayRate*timeNow);
	double freelivingNumber = freeliving*expo + (eggsProductionRate*(1.0-expo)/ReservoirDecayRate); // TODO: CHECK HERE

	return freelivingNumber;
}
