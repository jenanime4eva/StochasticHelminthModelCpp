/*
Individual-based stochastic simulation of helminths
Written by Jie Yang, November 2014
Compile this specific file to generate the .exe file
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <random>

int main(int argc, char** argv)
{
	// Do we have the correct command line?
	if(argc!=4)
	{
		std::cout << "Incorrect command line. Expected:\nFile.exe logFile paramFile resultsFile\n";
		return 0; 
	}


	return 1; 
}
