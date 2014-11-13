/*
Individual-based stochastic simulation of helminths
Written by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <random>

#include "CParamReader.h"

//using namespace std;

#define BUFFER_SIZE 1024

int main(int argc, char** argv)
{
	// Do we have the correct command line?
	/*
	if(argc!=2)
	{
		std::cout << "Incorrect command line. Expected:\nFile.exe logFile paramFile resultsFile\n";
		return 0; 
	}
	*/

	// Debugging stuff...
	// spit out the command line...
	for(int i=0;i<argc;i++)
	{
		std::cout << argv[i] << "\n";
	}

	CParamReader myReader;

	if(!myReader.setFileName(argv[1]))
	{
		std::cout << "Failed to attach file: " << argv[1] << "\n";
	} else
	{
		std::cout << "Attached file: " << argv[1] << "\n";
	}

	char* test = myReader.geteParamString("param1");
	test = myReader.geteParamString("param2");

	return 1; 
}
