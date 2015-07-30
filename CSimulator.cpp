/*
 * File Name: CSimulator.cpp
 *
 *  Created on: 26 Nov 2014
 *  Author: Jie Yang
 *
 */

#include "CSimulator.h"
#include "CHost.h"
#include "CParamReader.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <vector.h>
#include "randlib.h"

using namespace std;

// Class Constructor
CSimulator::CSimulator() {
	// Set default values

	// Model running parameters
	nRepetitions = 0;
	nYears = 0;
	nHosts = 0;

	// Social structure
	contactAgeBreaks = NULL;
	contactAgeBreaksLength = 0;
	betaValues = NULL;
	betaValuesLength = 0;
	rhoValues = NULL;
	rhoValuesLength = 0;

	// Demographic structure
	survivalCurve = survivalCurveCumul = NULL;
	hostMu = probDeath = probDeathIntegral = NULL;
	maxHostAgeCompare = NULL;
	demogDt = 0;
	hostMuData = NULL;
	hostMuDataLength = 0;
	muDataUpperBounds = NULL;
	muDataUpperBoundsLength = 0;
	upperAgeBound = 0;
	maxDtIntervals = 0;
	maxHostAge = 0;
	tinyIncrement = 0.01;

	// Epidemiological parameters
	k = 0;
	lambda = 0;
	R0 = 0;
	ReservoirDecayRate = 0;
	sigma = 0;
	gamma = 0;
	z = 0;
	psi = 0;

	// Treatment parameters
	treatmentBreaks = 0;
	treatmentBreaksLength = 0;
	coverage = NULL;
	coverageLength = 0;
	drugEff = 0;
	treatStart = 0;
	treatEnd = 0;
	treatInterval = 0;
	treatmentTimes = NULL;
	treatmentTimesLength = 0;

	// Results
	surveyTimesDt = 0;
	surveyResultTimes = NULL;
	surveyResultTimesLength = 0;
	meanFemaleWorms = 0;
	meanInfantFemaleWorms = 0;
	meanPreSACFemaleWorms = 0;
	meanSACFemaleWorms = 0;
	meanAdultFemaleWorms = 0;

	// Auxiliary functions
	vectorArray = NULL;
	contactIndices = NULL;
	treatmentIndices = NULL;

	myRealization = NULL;
}

// Class Destructor
CSimulator::~CSimulator() {
	if (contactAgeBreaks != NULL)
		delete[] contactAgeBreaks;

	if (betaValues != NULL)
		delete[] betaValues;

	if (rhoValues != NULL)
		delete[] rhoValues;

	if (hostMuData != NULL)
		delete[] hostMuData;

	if (muDataUpperBounds != NULL)
		delete[] muDataUpperBounds;

	if (maxHostAgeCompare != NULL)
		delete[] maxHostAgeCompare;

	if (survivalCurve != NULL)
		delete[] survivalCurve;

	if (survivalCurveCumul != NULL)
		delete[] survivalCurveCumul;

	if (hostMu != NULL)
		delete[] hostMu;

	if (probDeath != NULL)
		delete[] probDeath;

	if (probDeathIntegral != NULL)
		delete[] probDeathIntegral;

	if (treatmentBreaks != NULL)
		delete[] treatmentBreaks;

	if (coverage != NULL)
		delete[] coverage;

	if (treatmentTimes != NULL)
		delete[] treatmentTimes;

	if (surveyResultTimes != NULL)
		delete[] surveyResultTimes;

	if (vectorArray != NULL)
		delete[] vectorArray;

	if (contactIndices != NULL)
		delete[] contactIndices;

	if (treatmentIndices != NULL)
		delete[] treatmentIndices;

	int r;
	if(myRealization!=NULL)
	{
		for(r=0;r<nRepetitions;r++)
			delete myRealization[r];

		delete[] myRealization;
	}
}

// Initialise the input/output aspects of the simulator
bool CSimulator::initialiseIO(char* run, char* path, char* paramFilePath)
{
	// As string classes
	runName = run;
	thePath = path;

	// Setting up all the files needed
	std::string logFilePath = thePath + runName + ".log.txt";
	logStream.open(logFilePath.c_str());
	if (!logStream.is_open()) {
		std::cout << "Couldn't open the log file: " << logFilePath
				<< "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Log file opened.\n" << std::flush;

	myReader.setNewFileName(paramFilePath);
	if (!myReader.setNewFileName(paramFilePath)) {
		logStream << "Couldn't open the parameter file: " << paramFilePath
				<< "\nexiting\n" << std::flush;
		return false;
	}
	logStream << "Param reader configured.\n" << std::flush;

	// Stub name for all results files
	resultsStub = thePath + runName;
	logStream << "Stub name for all results files: " << resultsStub << "\n" << std::flush;

	// Everything is ok
	return true;
}

// General initialisation
bool CSimulator::initialiseSimulation()
{
	char* temp; // General purpose pointer to string

	// READ IN MODEL RUNNING PARAMETERS

	logStream << "\nTEST PARAMETERS READ IN\n\n" << std::flush;

	// Number of repetitions
	nRepetitions = atoi(myReader.getParamString("repNum"));
	logStream << "Number of repetitions: " << nRepetitions << "\n" << std::flush; // Test flag

	// Number of years to run
	nYears = atoi(myReader.getParamString("nYears"));
	logStream << "Number of years to run: " << nYears << "\n" << std::flush; // Test flag

	// Number of hosts
	nHosts = atoi(myReader.getParamString("nHosts"));
	logStream << "Number of hosts: " << nHosts << "\n" << std::flush; // Test flag

	// SET UP SOCIAL STRUCTURE

	// Read in contact age group breaks
	temp = myReader.getParamString("contactAgeBreaks");
	if (temp != NULL) {
		contactAgeBreaks = readDoublesVector(temp, contactAgeBreaksLength);
	}

	// Read in beta values (contact rates)
	temp = myReader.getParamString("betaValues");
	if (temp != NULL) {
		betaValues = readDoublesVector(temp, betaValuesLength);
	}

	// Read in rho values (contribution to the reservoir by contact age group)
	temp = myReader.getParamString("rhoValues");
	if (temp != NULL) {
		rhoValues = readDoublesVector(temp, rhoValuesLength);
	}

	/// SET UP DEMOGRAPHY

	// Read in host death rates
	temp = myReader.getParamString("hostMuData");
	if (temp != NULL) {
		hostMuData = readDoublesVector(temp, hostMuDataLength);
	}

	// Read in host death rate upper bounds
	temp = myReader.getParamString("upperBoundData");
	if (temp != NULL) {
		muDataUpperBounds = readDoublesVector(temp, muDataUpperBoundsLength);
	}

	// Age step for survival curve in years
	demogDt = atof(myReader.getParamString("demogDt"));

	// Construct mu, probability of death and survival vectors
	int currentMuIndex = 0;
	double currentSurvival = 1;
	double currentSurvivalCumul = 0;
	double currentProbDeathCumul = 0;
	double currentMuDtCumul = 0;

	maxDtIntervals = (int) floor(muDataUpperBounds[muDataUpperBoundsLength-1]/demogDt);
	upperAgeBound = maxDtIntervals * demogDt;

	maxHostAgeCompare = new double[2];
	maxHostAgeCompare[0] = upperAgeBound;
	maxHostAgeCompare[1] = contactAgeBreaks[contactAgeBreaksLength-1];
	maxHostAge = (int) min(maxHostAgeCompare,2); 	// Get maximum host age
	logStream << "\nmaxHostAge: " << maxHostAge << "\n";

	// Recalculate maxDtIntervals now that we know the maxHostAge value
	maxDtIntervals = maxHostAge/demogDt;

	survivalCurve = new double[maxDtIntervals];
	survivalCurveCumul = new double[maxDtIntervals];
	hostMu = new double[maxDtIntervals];
	probDeath = new double[maxDtIntervals];
	probDeathIntegral = new double[maxDtIntervals];

	for (int i=0;i<maxDtIntervals;i++)
	{
		double currentIntEnd = (i + 1) * demogDt;
		// Is current dt interval within data upper bound?
		if (muDataUpperBounds[currentMuIndex] + tinyIncrement < currentIntEnd)
		{
			currentMuIndex++; // Add one to currentMuIndex
		}
		hostMu[i] = hostMuData[currentMuIndex];

		probDeath[i] = currentSurvival * hostMu[i] * demogDt;
		probDeathIntegral[i] = currentProbDeathCumul + probDeath[i];
		currentProbDeathCumul = probDeathIntegral[i]; // Cumulative probability of dying in the ith year

		currentMuDtCumul += hostMu[i] * demogDt;
		survivalCurve[i] = exp(-currentMuDtCumul);
		currentSurvival = survivalCurve[i]; // Host survival curve

		survivalCurveCumul[i] = survivalCurve[i] + currentSurvivalCumul;
		currentSurvivalCumul = survivalCurveCumul[i];
	}

	// READ IN EPIDEMIOLOGICAL PARAMETERS

	// Shape parameter of assumed negative binomial distribution of worms amongst host
	k = atof(myReader.getParamString("k"));

	// Eggs per gram
	lambda = atof(myReader.getParamString("lambda"));

	// Basic reproductive number
	R0 = atof(myReader.getParamString("R0"));
	logStream << "\nR0: " << R0 << "\n" << std::flush; // Test flag

	// Decay rate of eggs in the environment
	ReservoirDecayRate = atof(myReader.getParamString("ReservoirDecayRate"));

	// Worm death rate i.e. 1/worm life span, same for all development stages
	sigma = atof(myReader.getParamString("sigma"));

	// Exponential density dependence of parasite adult stage (N.B. fecundity parameter z = exp(-gamma))
	gamma = atof(myReader.getParamString("gamma"));
	z = exp(-gamma); // Fecundity parameter

	// Dummy psi value prior to R0 calculation
	psi = calculatePsi();
	logStream << "\npsi: " << psi << "\n" << std::flush; // Test flag

	// SET UP TREATMENT

	// Read in treatment age group breaks
	temp = myReader.getParamString("treatmentBreaks");
	if (temp != NULL) {
		treatmentBreaks = readDoublesVector(temp, treatmentBreaksLength);
	}
	logStream << "\ntreatmentBreaks vector length: " << treatmentBreaksLength << "\n" << std::flush; // Test flag
	logStream << "Infant treatment break: " << treatmentBreaks[0] << "\n" << std::flush; // Test flag
	logStream << "pre-SAC treatment break: " << treatmentBreaks[1] << "\n" << std::flush; // Test flag
	logStream << "SAC treatment break: " << treatmentBreaks[2] << "\n" << std::flush; // Test flag
	logStream << "Adult treatment break: " << treatmentBreaks[3] << "\n" << std::flush; // Test flag
	logStream << "Maximum treatment age: " << treatmentBreaks[4] << "\n" << std::flush; // Test flag

	// Read in coverages
	temp = myReader.getParamString("coverage");
	if (temp != NULL) {
		coverage = readDoublesVector(temp, coverageLength);
	}
	logStream << "\nCoverages vector length: " << coverageLength << "\n" << std::flush; // Test flag
	logStream << "Infant coverage: " << coverage[0] << "\n" << std::flush; // Test flag
	logStream << "pre-SAC coverage: " << coverage[1] << "\n" << std::flush; // Test flag
	logStream << "SAC coverage: " << coverage[2] << "\n" << std::flush; // Test flag
	logStream << "Adult coverage: " << coverage[3] << "\n\n" << std::flush; // Test flag

	// Drug efficacy
	drugEff = atof(myReader.getParamString("drugEff"));

	// Treatment year start
	treatStart = atoi(myReader.getParamString("treatStart"));
	// Treatment year end
	treatEnd = atoi(myReader.getParamString("treatEnd"));
	// Interval between treatments in years
	treatInterval = atof(myReader.getParamString("treatInterval"));
	// Now create a vector of treatment times
	treatmentTimesLength = (int) ((treatEnd - treatStart)/treatInterval)+1;
	treatmentTimes = new double[treatmentTimesLength];
	treatmentTimes[0] = treatStart; // First vector entry
	for (int i=1;i<treatmentTimesLength;i++)
	{
		treatmentTimes[i] = treatmentTimes[i-1] + treatInterval;
	}

	// SET UP RESULTS COLLECTION

	// Time step for the survey times in years
	surveyTimesDt = atof(myReader.getParamString("surveyTimesDt"));
	// Now create a vector of survey result times
	surveyResultTimesLength = (int) (nYears/surveyTimesDt)+1;
	surveyResultTimes = new double[surveyResultTimesLength];
	surveyResultTimes[0] = 0; // First vector entry
	for (int i=1;i<surveyResultTimesLength;i++)
	{
		surveyResultTimes[i] = surveyResultTimes[i-1] + surveyTimesDt;
	}

	// Set up arrays to be used later
	contactIndices = new int[maxHostAge];
	treatmentIndices = new int[maxHostAge];

	// Set up realisation array
	myRealization = new CRealization*[nRepetitions];

	// Set up realisations
	for (int repNo=0;repNo<nRepetitions;repNo++)
	{
		myRealization[repNo] = new CRealization;
		myRealization[repNo]->initialize(this,repNo);
	}

	return true;
}

// Run simulation
void CSimulator::runSimulation()
{
	// Loop through the runs
	for (int repNo=0;repNo<nRepetitions;repNo++)
	{
		// Run a realisation
		myRealization[repNo]->run(repNo);
	}
}

// Output simulation
void CSimulator::outputSimulation()
{
	// Output survey-type results if there's any
	if (surveyResultTimes != NULL)
	{
		// There's some survey results
		// Create output stream
		std::string surveyResultsOut = thePath + runName + ".surveyResults.txt";
		std::ofstream surveyStream(surveyResultsOut.c_str());

		// For printing values of the means of the four age groups over time
		for (int j = 0; j < surveyResultTimesLength; j++)
		{
			// Print times in the first column
			surveyStream << surveyResultTimes[j] << "\t";

			double sumFemaleWorms = 0;
			double sumInfantFemaleWorms = 0;
			double sumPreSACFemaleWorms = 0;
			double sumSACFemaleWorms = 0;
			double sumAdultFemaleWorms = 0;

			for (int repNo = 0; repNo < nRepetitions; repNo++)
			{
				// Do some calculations
				sumInfantFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][repNo].meanInfantFemaleWormsPerRun;
				sumPreSACFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][repNo].meanPreSACFemaleWormsPerRun;
				sumSACFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][repNo].meanSACFemaleWormsPerRun;
				sumAdultFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][repNo].meanAdultFemaleWormsPerRun;
				sumFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][repNo].meanFemaleWormsPerRun;
			}

			// Mean of the realisations
			meanInfantFemaleWorms = sumInfantFemaleWorms/nRepetitions;
			meanPreSACFemaleWorms = sumPreSACFemaleWorms/nRepetitions;
			meanSACFemaleWorms = sumSACFemaleWorms/nRepetitions;
			meanAdultFemaleWorms = sumAdultFemaleWorms/nRepetitions;
			meanFemaleWorms = sumFemaleWorms/nRepetitions;

			// Print out infant, pre-SAC, SAC, adult and whole population mean female worm burdens over survey time intervals
			surveyStream << meanInfantFemaleWorms << "\t" << std::flush;
			surveyStream << meanPreSACFemaleWorms << "\t" << std::flush;
			surveyStream << meanSACFemaleWorms << "\t" << std::flush;
			surveyStream << meanAdultFemaleWorms << "\t" << std::flush;
			surveyStream << "\t" << meanFemaleWorms << "\t" << std::flush;

			surveyStream << "\n" << std::flush;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
/// Auxiliary function definitions

// Creates a vector of doubles from string from the parameter file
double* CSimulator::readDoublesVector(char* currentString, int& currentVectorLength)
{
	// IMPORTANT: Make sure the string read in does not have trailing zeros
	// These are currently removed in the parameter reading function

	char* endPointer; // endPointer for each call to strto_ type functions

	// This part is just to count the number of entries in the vector
	int counter = 1;
	double* tempVector = new double[strlen(currentString)]; // Create temporary vector
	tempVector[0] = strtod(currentString, &endPointer);
	while(strlen(endPointer)>0)
	{
		tempVector[counter] = strtod(endPointer, &endPointer);
		counter++; // Add one to counter
	}
	// Delete tempVector as we don't need it anymore
	if (tempVector!=NULL)
		delete[] tempVector;

	// NOW create the vector of doubles of values we actually require
	vectorArray = new double[counter];
	vectorArray[0] = strtod(currentString, &endPointer);
	for(int i=1;i<counter;i++)
	{
		vectorArray[i] = strtod(endPointer, &endPointer);
	}

	currentVectorLength = counter; // Count number of elements in the vector

	return vectorArray;
}

// This function takes a random number (0-1) and multiplies it by the MAX of the passed array.
// It then finds the smallest index that has an array value greater than the product above.
// For a cumulative multinomial array, this will return the index of the event that occurred.
// Also used for drawing a lifespan from the survival curve integral.
int CSimulator::multiNomBasic1(double* array, int length, double randNum)
{
	int loopMax = ceil(log(length) / log(2) + 2);
	int bottom = -1;
	//double bottomVal;
	int top = length - 1;
	double topVal = array[top];
	double target = topVal * randNum;
	int count = 0;

	while (++count < loopMax && (top - bottom > 1)) {
		int mid = (top + bottom) / 2;
		double midVal = array[mid];
		if (midVal >= target) {
			top = mid;
			topVal = midVal;
		} else {
			bottom = mid;
			//bottomVal = midVal;
		}
	}

	if (count >= loopMax) {
		logStream << "Max iterations exceeded in multiNomBasic1(...),\n"
				<< std::flush;
		return -1;
	}

	return top;
}

// This function takes a random number (0-1) and multiplies it by the SUM of the passed array.
// It then finds the smallest index that has an array value greater than the product above.
// For a cumulative multinomial array, this will return the index of the event that occurred.
int CSimulator::multiNomBasic2(double* array, int length, double randNum)
{
	int loopMax = ceil(log(length) / log(2) + 2);
	int bottom = -1;
	int top = length-1;
	double topVal = sumArray(array,length);
	double target = topVal * randNum;
	int count = 0;

	while (++count < loopMax && (top - bottom > 1)) {
		int mid = (top + bottom) / 2;
		double midVal = sumArray(array,mid+1);
		if (midVal >= target) {
			top = mid;
			topVal = midVal;
		} else {
			bottom = mid;
		}
	}

	if (count >= loopMax) {
		logStream << "Max iterations exceeded in multiNomBasic2(...),\n"
				<< std::flush;
		return -1;
	}

	return top;
}

// Calculate the psi value
double CSimulator::calculatePsi()
{
	// Higher resolution
	double deltaT = 0.1;

	// hostMu for the new age intervals
	int currentHostMuGroupIndex = 0;
	int currentModelAgeGroupCatIndex = 0;
	double currentMeanDeathsCumul = 0;
	double sumHostSurvival = 0;
	double meanLifespan = 0;
	double sumRhoAgeHostSurvivalCurveK = 0;

	int maxHostAgeInterval = (int) (maxHostAge/deltaT)+1;
	double* modelAges = new double[maxHostAgeInterval];
	double* hostMuArray = new double[maxHostAgeInterval];
	double* hostSurvivalCurve = new double[maxHostAgeInterval];
	double* survivalCurveSum = new double[maxHostAgeInterval];
	double* betaAge = new double[maxHostAgeInterval];
	double* rhoAge = new double[maxHostAgeInterval];
	double* wSurvival = new double[maxHostAgeInterval];
	double* K = new double[maxHostAgeInterval];

	double ADD = 0;
	for (int i=0;i<maxHostAgeInterval;i++)
	{
		// Interval-centered ages for the age intervals
		modelAges[i] = min(contactAgeBreaks,sizeof(contactAgeBreaks)/sizeof(contactAgeBreaks[0])) + 0.5*deltaT + ADD;
		ADD = ADD + deltaT;

		double currentIntEnd = (i + 1)*deltaT*demogDt;
		if (muDataUpperBounds[currentHostMuGroupIndex] + tinyIncrement < currentIntEnd)
		{
			currentHostMuGroupIndex++; // Add one to currentHostMuGroupIndex
		}
		hostMuArray[i] = hostMuData[currentHostMuGroupIndex];

		currentMeanDeathsCumul += hostMuArray[i]*deltaT*demogDt;
		hostSurvivalCurve[i] = exp(-currentMeanDeathsCumul);

		survivalCurveSum[i] = hostSurvivalCurve[i] + sumHostSurvival;
		sumHostSurvival = survivalCurveSum[i];

		// Need rho and beta at this age resolution as well
		if (contactAgeBreaks[currentModelAgeGroupCatIndex] + tinyIncrement < currentIntEnd)
		{
			currentModelAgeGroupCatIndex++; // Add one to currentModelAgeGroupCatIndex
		}
		betaAge[i] = betaValues[currentModelAgeGroupCatIndex-1];
		rhoAge[i] = rhoValues[currentModelAgeGroupCatIndex-1];

		wSurvival[i] = exp(-sigma*modelAges[i]);
	}

	meanLifespan = sumHostSurvival*deltaT;

	K[0] = betaAge[0]*wSurvival[0]*deltaT; // First K array entry
	for(int n=1;n<maxHostAgeInterval;n++)
	{
		double sumBetaAgewSurvival = 0;
		int a = n;
		for(int j=0;j<=n;j++)
		{
			sumBetaAgewSurvival += betaAge[j]*wSurvival[a];
			a = a-1;
		}
		K[n] = sumBetaAgewSurvival*deltaT;
		sumRhoAgeHostSurvivalCurveK += rhoAge[n]*hostSurvivalCurve[n]*K[n];
	}

	double summation = sumRhoAgeHostSurvivalCurveK*deltaT;

	psi = R0*meanLifespan*ReservoirDecayRate/(lambda*z*summation);

	// Delete arrays we don't need anymore to free up memory
	if (modelAges!=NULL)
		delete[] modelAges;

	if (hostMuArray!=NULL)
		delete[] hostMuArray;

	if (hostSurvivalCurve!=NULL)
		delete[] hostSurvivalCurve;

	if (survivalCurveSum!=NULL)
		delete[] survivalCurveSum;

	if (betaAge!=NULL)
		delete[] betaAge;

	if (rhoAge!=NULL)
		delete[] rhoAge;

	if (wSurvival!=NULL)
		delete[] wSurvival;

	if (K!=NULL)
		delete[] K;

	return (2*psi); // This 2 is here because it's infection with both male and female worms and only half of these are going to be female
}

// Draw a life span from the survival curve from the population
double CSimulator::drawLifespan()
{
	// Get a random integer from the probDeathIntegral using the multinomial generator. This shouldn't be zero!
	double currentRand = myRandUniform();

	int index = multiNomBasic1(probDeathIntegral, maxDtIntervals, currentRand); // maxDtIntervals = maxHostAge here?

	// Choose a point in the middle of the interval in which the person dies
	double ans = (index + 0.5) * demogDt;

	return ans;
}

// Contact index array
int* CSimulator::contactAgeGroupIndex()
{
	int currentContactIndex = 1;
	for(int i=0;i<maxHostAge;i++)
	{
		contactIndices[i] = currentContactIndex-1;
		if(i >= contactAgeBreaks[currentContactIndex])
		{
			currentContactIndex++;
		}
	}

	return contactIndices;
}

// Treatment level index array
int* CSimulator::treatmentAgeGroupIndex()
{
	int currentTreatIndex = 1;
	for(int i=0;i<maxHostAge;i++)
	{
		treatmentIndices[i] = currentTreatIndex-1;
		if(i >= treatmentBreaks[currentTreatIndex])
		{
			currentTreatIndex++;
		}
	}

	return treatmentIndices;
}

// Uniform distribution random number generator (generates real number between 0 to 1)
double CSimulator::myRandUniform()
{
	//genunf(double low, double high) generates uniform real between low and high (see randlib files)
	  return  genunf(0,1);
}

// Gamma distribution random number generator (for returning the value of si, see myRandPoisson function)
double CSimulator::myRandGamma(double l, double s)
{
	// double gengam(double l,double s) generates random deviates from gamma distribution (see randlib files)
	// l is location parameter, s is the shape parameter
	return gengam(l,s);
}

// Poisson distribution random number generator (for returning total worm numbers)
double CSimulator::myRandPoisson(double mu)
{
    // long ignpoi(double mu) generates a single random deviate from a poisson distribution with mean mu (see randlib files)
	return ignpoi(mu);
}

// Binomial distribution random number generator
double CSimulator::myRandBinomial(long n, double p)
{
    // long ignbin(long n,double p) generates a single random deviate from a binomial distribution (see randlib files)
	// n is the number of trials, p is probability of event
	return ignbin(n,p);
}

// Exponential distribution random number generator
double CSimulator::myRandExponential(double rateValue)
{
    // double genexp(double a) generates a single random deviate from an exponential distribution with mean a (see randlib files)
	// N.B. Mean = 1/rateValue
	return genexp(1/rateValue);
}

// Find the sum of an array
double CSimulator::sumArray(double* array,int arrayLength)
{
	double sum = 0;
	for(int n=0;n<arrayLength;n++)
	{
	  sum += array[n];
	}

	return sum;
}

// Function to find minimum of a list of values
double CSimulator::min(double* Numbers, int Count)
{
	double Minimum = Numbers[0];

	for(int i=0;i<Count;i++)
	{
		if(Minimum>Numbers[i])
		{
			Minimum = Numbers[i];
		}
	}

	return Minimum;
}

// Function to find maximum of a list of values
double CSimulator::max(double* Numbers, int Count)
{
	double Maximum = Numbers[0];

	for(int i=0;i<Count;i++)
	{
		if(Maximum<Numbers[i])
		{
			Maximum = Numbers[i];
		}
	}

	return Maximum;
}

// Return index of smallest element in an array
int CSimulator::indexSmallestElement(double* array, int size)
{
    int smallestIndex = 0;
    double temp=array[0];
    int i;
    for(i=0;i<size;i++)
    {
        if(array[i]<temp)
        {
            smallestIndex = i;
            temp=array[i];
        }
    }

    return smallestIndex;
}

