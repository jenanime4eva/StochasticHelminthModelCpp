/*
 * CParamReader.cpp
 *
 *  Created on: 12 Nov 2014
 *      Author: jtrusc
 */

#include "CParamReader.h"
#include <string.h>

CParamReader::CParamReader()
{
	paramBuffer = new char[BUFFER_SIZE];
	filePathString = NULL;
}

CParamReader::~CParamReader()
{
	delete[] paramBuffer;
	if(filePathString!=NULL)
		delete[] filePathString;
}


// set the file name and check that it exists.
bool CParamReader::setFileName(char* filePath)
{
	paramFileStream.open(filePath);

	if(!paramFileStream.is_open())
		return false;

	paramFileStream.close(); // just checking; close it.

	//filePathString = new char[strlen(filePath)+1];
	//strcpy(filePathString, filePath);
	filePathString = filePath;

	return true;
}

// return a string containing parameter data or Null.
char* CParamReader::geteParamString(const char* paramName)
{
	// do we have a file attached?
	paramFileStream.open(filePathString);

	if(!paramFileStream.is_open())
		return false;

	int paramLength = strlen(paramName);

	bool found = false;
	char* paramString = NULL;
	char* token;
	// run through the lines of the file for the parameter.
	while(!paramFileStream.eof() && !found)
	{
		paramFileStream.getline(paramBuffer,BUFFER_SIZE-1); // I don't know about the final \0.
		token = strtok(paramBuffer,"\t");
		// anything on this line?
		if(token!=NULL)
		{
			// line doesn't begin with # and matches param name?
			if(token[0]!='#' && strncmp(paramName,token,paramLength)==0)
			{
				// found it!
				found = true;
				paramString = strtok(NULL,"\t");
			}
		}
	}

	paramFileStream.close();

	return paramString;
}



