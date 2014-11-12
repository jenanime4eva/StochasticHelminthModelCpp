/*
 * CParamReader.h
 *
 *  Created on: 12 Nov 2014
 *      Author: jtrusc
 */

#ifndef CPARAMREADER_H_
#define CPARAMREADER_H_

#include <fstream>

#define BUFFER_SIZE 1024

class CParamReader
{
public:
	CParamReader();
	virtual ~CParamReader();

	bool setFileName(char* filePath);
	char* geteParamString(char* paraName);


	// members...
	std::ifstream paramFileStream;
	char* paramBuffer;
	char* filePathString;

	// use getline for the line getting.
};

#endif /* CPARAMREADER_H_ */
