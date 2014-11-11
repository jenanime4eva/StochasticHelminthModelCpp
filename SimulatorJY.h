#pragma once //  Causes the current source file to be included only once in a single compilation, may improve compilation speed

#include <fstream>

#define BUFFER_SIZE 1024

struct meanWormBurden
{
	double time;
	int M;
};

class helsimJY
{
public:
	helsimJY(void); // Declare constructor…
	virtual ~helsimJY(void); //…and destructor

	// File input and output…
	std::ifstream paramsIn; // Input file stream for parameters
	std::ofstream logStream, resultsStream; // Output to log file and results
	bool initialiseIO(char* logFileName, char* paramFileName, char* resultsFileName); 
private:
	bool getNextLine(std::ifstream* currentStream); // Get next usable line
	bool isLineRubbish(char* currentLine); 
public:
	// Initialisation...
	char buffer[BUFFER_SIZE]; // A buffer is a region of a physical memory storage used to temporarily store data while it is being moved from one place to another
	bool initialiseSimulation(); 

	// Simulation...
	// Host stuff (double type values)
	double demogEta, demogB, InfantBeta, PreSACBeta, SACBeta, AdultBeta;
	// Worm stuff (double type values)
	double R0, gamma, z, k, sigma, LDecayRate;
	// Treatment stuff (double type values)
	double InfantCoverage, PreSACCoverage, SACCoverage, AdultCoverage, drugEfficacy, treatFreq;
	// Simulation stuff (double type values)
	double tstep, timeElapsed;
	// Host stuff (integer type values)
	int hostPopulation, nAgeGroups, contactBreakInfants, contactBreakPreSAC, contactBreakSAC, contactBreakAdults;
	// Worm stuff (integer type values)
	int lambda;
	// Treatment stuff (integer type values)
	int treatBreakInfants, treatBreakPreSAC, treatBreakSAC, treatBreakAdults, treatStartYear, treatEndYear;
	// Simulation stuff (integer type values)
	int nRepetitions, nMaxYear, nOutFreq, nTimeSteps;
	meanWormBurden** results; // Array of results
	void runSimulation();

	// Results output and processing...
	void outputSimulation(int n);
	void outputWormMean();
	void outputWormVariance();
};
