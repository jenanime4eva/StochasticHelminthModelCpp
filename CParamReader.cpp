/*
 * CParamReader.cpp
 *
 *  Created on: 12 Nov 2014
 *      Author: jtrusc
 *  Modified: 26 Nov 2014
 *		Author: Jie Yang
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


// Set the file name and check that it exists
bool CParamReader::setNewFileName(char* filePath)
{
	// If there's a file attached, remove it.
	if(paramFileStream.is_open())
		paramFileStream.close();

	paramFileStream.open(filePath);

	if(!paramFileStream.is_open())
		return false;

	paramFileStream.close(); // Just checking, close it

	//filePathString = new char[strlen(filePath)+1];
	//strcpy(filePathString, filePath);
	filePathString = filePath;

	return true;
}

// Return a string containing parameter data or Null
char* CParamReader::getParamString(const char* paramName)
{
	// Do we have a file attached?
	paramFileStream.open(filePathString);

	if(!paramFileStream.is_open())
		return NULL;

	int paramLength = strlen(paramName);

	bool found = false;
	char* paramString = NULL;
	char* token;
	// Run through the lines of the file for the parameter
	while(!paramFileStream.eof() && !found)
	{
		paramFileStream.getline(paramBuffer,BUFFER_SIZE-1);
		token = strtok(paramBuffer,"\t");
		// Anything on this line?
		if(token!=NULL)
		{
		// Line doesn't begin with # and matches param name?
			if(token[0]!='#' && strncmp(paramName,token,paramLength)==0)
			{
			// Found it!
			found = true;
			paramString = strtok(NULL,"\t");
			}
		}
	}

	paramFileStream.close();

	return paramString;
}



