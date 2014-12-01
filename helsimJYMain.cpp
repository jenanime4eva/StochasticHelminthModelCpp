/*
Individual-based stochastic simulation of helminths
Written by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/


//#include "CParamReader.h"
#include "CSimulator.h"
//#include "CHost.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <random>

// Using namespace std;

#define BUFFER_SIZE 1024

int main(int argc, char** argv)
{
	// Do we have the correct command line?

	if(argc!=4)
	{
		std::cout << "Incorrect command line. Expected arguments for your .exe file: logFile paramFile resultsFile\n";
		return 0; 
	}

	CSimulator Simulate;
	Simulate.initialiseIO(argv[1],argv[2],argv[3]);

	// Debugging stuff...
	// Spit out the command line...
//	for(int i=0;i<argc;i++)
//	{
//		std::cout << argv[i] << "\n";
//	}
//
//	CParamReader myReader;
//
//	if(!myReader.setNewFileName(argv[1]))
//	{
//		std::cout << "Failed to attach file: " << argv[1] << "\n";
//	} else
//	{
//		std::cout << "Attached file: " << argv[1] << "\n";
//	}

	//char* test = myReader.getParamString("param1");
	//test = myReader.getParamString("param2");

	return 1; 
}
