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
	treatmentOnOff = 0;
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

	if(myRealization!=NULL)
	{
		for(int r=0;r<nRepetitions;r++)
			delete myRealization[r];

		delete[] myRealization;
	}
}

// Initialise the input/output aspects of the simulator
bool CSimulator::initialiseIO(char* run, char* path, char* paramFilePath)
{
	// As string classes

	printf("\n"); // I don't know why but removing this printf statement causes the code to crash (runs fine in debug mode though)
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
	nYears = atof(myReader.getParamString("nYears"));
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

	// Is treatment on or off?
	treatmentOnOff = atoi(myReader.getParamString("treatmentOnOff"));

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
	drugEff = (double) drugEff*treatmentOnOff;

	// Treatment year start
	treatStart = atof(myReader.getParamString("treatStart"));
	// Treatment year end
	treatEnd = atof(myReader.getParamString("treatEnd"));
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

	// Set up arrays to be used later
	contactIndices = new int[maxHostAge+1];
	treatmentIndices = new int[maxHostAge+1];

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
bool CSimulator::runSimulation()
{
	// Loop through the runs
	for (int repNo=0;repNo<nRepetitions;repNo++)
	{
		// Run a realisation
		myRealization[repNo]->run(repNo);
	}

	return true;
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

		// Uncomment as necessary for what kind of output file you want

		///* TODO: FOR TIME COURSE PLOTS OF FEMALE WORM BURDEN OF THE INFANT, PRE-SAC, SAC AND ADULT AGE GROUPS
		for (int j = 0; j < surveyResultTimesLength; j++)
		{
			// Print times in the first column
			surveyStream << surveyResultTimes[j] << "\t";

			// These variables reset after each iteration
			int sumInfantFemaleWorms = 0;
			int sumPreSACFemaleWorms = 0;
			int sumSACFemaleWorms = 0;
			int sumAdultFemaleWorms = 0;
			int sumFemaleWorms = 0;
			int sumInfantNumber = 0;
			int sumPreSACNumber = 0;
			int sumSACNumber = 0;
			int sumAdultNumber = 0;
			int sumTotalHostNumber = 0;

			// Loop through the hosts...
			for(int i=0;i<nHosts;i++)
			{
				// Loop through the runs...
				for(int repNo=0;repNo<nRepetitions;repNo++)
				{
					// Perform some calculations here

					// Infants
					sumInfantFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][i].infantFemaleWorms;
					sumInfantNumber += myRealization[repNo]->surveyResultsArrayPerRun[j][i].infantNumber;

					// Pre-SAC
					sumPreSACFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][i].preSACFemaleWorms;
					sumPreSACNumber += myRealization[repNo]->surveyResultsArrayPerRun[j][i].preSACNumber;

					// SAC
					sumSACFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][i].SACFemaleWorms;
					sumSACNumber += myRealization[repNo]->surveyResultsArrayPerRun[j][i].SACNumber;

					// Adults
					sumAdultFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][i].adultFemaleWorms;
					sumAdultNumber += myRealization[repNo]->surveyResultsArrayPerRun[j][i].adultNumber;

					// Whole population
					sumFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][i].femaleWorms;
					sumTotalHostNumber += myRealization[repNo]->surveyResultsArrayPerRun[j][i].totalHostNumber;
				}
			}

			// Print out infant, pre-SAC, SAC, adult and whole population mean female worm burdens over time

			// Infants
			surveyStream << (double) sumInfantFemaleWorms/(sumInfantNumber+0.01) << "\t"; // 0.01 is to avoid division by zero

			// Pre-SAC
			surveyStream << (double) sumPreSACFemaleWorms/(sumPreSACNumber+0.01) << "\t"; // 0.01 is to avoid division by zero

			// SAC
			surveyStream << (double) sumSACFemaleWorms/(sumSACNumber+0.01) << "\t"; // 0.01 is to avoid division by zero

			// Adults
			surveyStream << (double) sumAdultFemaleWorms/(sumAdultNumber+0.01) << "\t"; // 0.01 is to avoid division by zero

			// Whole Population
			surveyStream << "\t" << (double) sumFemaleWorms/(sumTotalHostNumber+0.01) << "\t"; // 0.01 is to avoid division by zero

			surveyStream << "\n"; // New line before next time output
		}
		//*/

		/* TODO: FOR TIME COURSE PLOTS OF FEMALE WORM BURDEN FOR EACH HOST (Set repNum = 1 in the parameter file)
		for (int j = 0; j < surveyResultTimesLength; j++)
		{
			// Print times in the first column
			surveyStream << surveyResultTimes[j] << "\t";

			// Loop through the hosts...
			for(int i=0;i<nHosts;i++)
			{
				// Loop through the runs
				for(int repNo=0;repNo<nRepetitions;repNo++)
				{
					// Print female worms for current host
					surveyStream << myRealization[repNo]->surveyResultsArrayPerRun[j][i].femaleWorms << "\t";
				}
			}

			surveyStream << "\n"; // New line before next time output
		}
		*/

		/* TODO: FOR WORKING OUT PROPORTION OF FEMALE WORM EXTINCTIONS VS YEARS AFTER TREATMENT
		int treatStartIndex = (int) (treatStart/surveyTimesDt);
		for (int j = treatStartIndex; j < surveyResultTimesLength; j++)
		{
			// Print times from first treatment time in the first column
			// N.B. Printout will actually be the current run time minus the treatment start year
			double timeOut = surveyResultTimes[j]-treatStart;
			surveyStream << timeOut << "\t";

			int sumZeroFemaleWorms = 0;

			// Loop through the hosts...
			for(int repNo=0;repNo<nRepetitions;repNo++)
			{
				int zeroFemaleWorms = 0;
				int sumFemaleWorms = 0;

				// Loop through the runs...
				for(int i=0;i<nHosts;i++)
				{
					sumFemaleWorms += myRealization[repNo]->surveyResultsArrayPerRun[j][i].femaleWorms;
				}

				// Has female worms reached zero yet? Allocate value of one for yes, zero for no
				if(sumFemaleWorms==0)
				{
					zeroFemaleWorms = 1; // Worms have gone extinct for current run
				}
				else
				{
					zeroFemaleWorms = 0; // Worms have not gone extinct for current run
				}

				sumZeroFemaleWorms += zeroFemaleWorms; // Sum up the number of runs that have reached zero female worms
			}

			// Output proportion of female worm extinctions out of the specified number of runs in the second column
			surveyStream << (double) sumZeroFemaleWorms/nRepetitions << "\n";
		}
		*/

		/* TODO: FOR WORM BURDEN VS AGE AT EQUILIBRIUM
		int equiTimeIndex = (int) (treatStart/surveyTimesDt)-1; // Get index of equilibrium time point
		for (int ageCat = 1; ageCat <= maxHostAge; ageCat++)
		{
			// Print host ages in the first column
			surveyStream << ageCat << "\t";

			// Reset these after each iteration
			int femaleWormsAgeCat = 0;
			int hostAge = 0;
			int hostCount = 0;

			// Loop through the hosts...
			for(int i=0;i<nHosts;i++)
			{
				// Loop through the runs...
				for(int repNo=0;repNo<nRepetitions;repNo++)
				{
					// Perform some calculations here

					// Age of host falls into current age category?
					hostAge = (int) ceil(myRealization[repNo]->surveyResultsArrayPerRun[equiTimeIndex][i].hostAge);
					if (hostAge == ageCat)
					{
						// Female worms at last survey time point
						femaleWormsAgeCat += myRealization[repNo]->surveyResultsArrayPerRun[equiTimeIndex][i].femaleWorms;

						// Add one to host count
						hostCount++;
					}
				}

			}

			surveyStream << (double) femaleWormsAgeCat/(hostCount+0.01); // 0.01 is to avoid division by zero

			surveyStream << "\n"; // New line before next age output
		}
		*/
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
int CSimulator::multiNomBasic(double* array, int length, double randNum)
{
	int loopMax = ceil(log(length) / log(2) + 2); // N.B. LOGb(x)/LOGb(a)=LOGa(x)
	int bottom = -1;
	int top = length - 1;
	double topVal = array[top];
	double target = topVal * randNum;
	int count = 0;

	while (++count < loopMax && (top - bottom > 1))
	{
		int mid = (top + bottom) / 2;
		double midVal = array[mid];

		if (midVal > target)
		{
			top = mid;
		} else
		{
			bottom = mid;
		}
	}

	if (count >= loopMax)
	{
		logStream << "Max iterations exceeded in multiNomBasic(...),\n"
				<< std::flush;
		return -1;
	}

	return top;
}

int CSimulator::multiNomBasicAlt(double* array, int length, double randNum)
{
	int loopMax = ceil(log(length) / log(2) + 2); // N.B. LOGb(x)/LOGb(a)=LOGa(x)
	int bottom = -1;
	int top = length - 1;
	double topVal = sumArray(array,top+1);
	double target = topVal * randNum;
	int count = 0;

	while (++count < loopMax && (top - bottom > 1))
	{
		int mid = (top + bottom) / 2;
		double midVal = sumArray(array,mid+1);

		if (midVal > target)
		{
			top = mid;
		} else
		{
			bottom = mid;
		}
	}

	if (count >= loopMax)
	{
		logStream << "Max iterations exceeded in multiNomBasicAlt(...),\n"
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

	int maxHostAgeInterval = (int) (maxHostAge/deltaT);
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
	double randNum = myRandUniform();
	int index = multiNomBasic(probDeathIntegral,maxDtIntervals,randNum); // maxDtIntervals = maxHostAge here?

	// Choose a point in the middle of the interval in which the person dies
	double ans = (index + 0.5) * demogDt;

	return ans;
}

// Contact index array
int* CSimulator::contactAgeGroupIndex()
{
	int currentContactIndex = 1;
	for(int i=0;i<=maxHostAge;i++)
	{
		if(i >= contactAgeBreaks[currentContactIndex])
		{
			currentContactIndex++;
		}

		contactIndices[i] = currentContactIndex-1;
	}

	return contactIndices;
}

// Treatment level index array
int* CSimulator::treatmentAgeGroupIndex()
{
	int currentTreatIndex = 1;
	for(int i=0;i<=maxHostAge;i++)
	{
		if(i >= treatmentBreaks[currentTreatIndex])
		{
			currentTreatIndex++;
		}

		treatmentIndices[i] = currentTreatIndex-1;
	}

	return treatmentIndices;
}

// Uniform distribution random number generator (generates real number between 0 to 1)
double CSimulator::myRandUniform()
{
	//genunf(double low, double high) generates uniform real between low and high (see randlib files)
	  return  genunf(0,1);
}

// Gamma distribution random number generator
double CSimulator::myRandGamma(double l, double s)
{
	// double gengam(double l,double s) generates random deviates from gamma distribution (see randlib files)
	// l is location parameter, s is the shape parameter
	return gengam(l,s);
}

// Poisson distribution random number generator
int CSimulator::myRandPoisson(double mu)
{
    // long ignpoi(double mu) generates a single random deviate from a poisson distribution with mean mu (see randlib files)
	return (int) ignpoi(mu);
}

// Binomial distribution random number generator
int CSimulator::myRandBinomial(long n, double p)
{
    // long ignbin(long n,double p) generates a single random deviate from a binomial distribution (see randlib files)
	// n is the number of trials, p is probability of event
	return (int) ignbin(n,p);
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
