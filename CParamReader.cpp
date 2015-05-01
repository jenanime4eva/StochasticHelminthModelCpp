/*
 * File Name: CParamReader.cpp
 *
 * Created on: 12 Nov 2014
 *
 */

#include "CParamReader.h"
#include <string.h>

// Class Constructor
CParamReader::CParamReader()
{
	paramBuffer = new char[BUFFER_SIZE];
	filePathString = NULL;
}

// Class Destructor
CParamReader::~CParamReader()
{
	delete[] paramBuffer;
	if(filePathString!=NULL)
		delete[] filePathString;
}

// Set up the file name and check that it exists
bool CParamReader::setNewFileName(char* filePath)
{
	// If there's a file attached, remove it.
	if(paramFileStream.is_open())
		paramFileStream.close();

	paramFileStream.open(filePath);

	if(!paramFileStream.is_open())
		return false;

	paramFileStream.close(); // Just checking, close it

	filePathString = new char[strlen(filePath)+1];
	strcpy(filePathString, filePath);

	return true;
}

// Return a string containing parameter data or NULL
char* CParamReader::getParamString(const char* paramName)
{
	// Do we have a file attached?
	paramFileStream.open(filePathString);

	if(!paramFileStream.is_open())
		return NULL;

	bool found = false;
	char* paramString = NULL;
	char* token;

	// Run through the lines of the file for the parameter
	while(!paramFileStream.eof() && !found) // eof = end of file function
	{
		paramFileStream.getline(paramBuffer,BUFFER_SIZE-1);
		token = strtok(paramBuffer,"\t");

		// Anything on this line?
		if(token!=NULL)
		{
			// Line doesn't begin with # and matches parameter name
			if(token[0]!='#' && strcmp(paramName,token)==0)
			{
				// Found it!
				found = true;
				paramString = strtok(NULL,"\t");
			}
		}
	}

	// Remove any extra spaces from the end of the string (as long as the string isn't zero length)
	int nullTerminal = strlen(paramString);
	while(paramString[nullTerminal-1]==' ')
		nullTerminal--;

	paramString[nullTerminal] = '\0'; // "\0" is the NULL character (it has the value 0)

	paramFileStream.close();

	return paramString;
}
