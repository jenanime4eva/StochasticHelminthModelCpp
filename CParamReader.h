/*
 * File name: CParamReader.h
 *
 * Created on: 12 Nov 2014
 *
 *      Each informative line of the parameter file should have the form:
 *      Name (TAB) Values separated by spaces (TAB) ## Comments
 *      Any line starting with '#' is completely ignored.
 *
 */

#ifndef CPARAMREADER_H_
#define CPARAMREADER_H_

#include <fstream>

#define BUFFER_SIZE 1024 // Define buffer size

class CParamReader
{
public:
	CParamReader(); // This is the constructor declaration
	virtual ~CParamReader(); // This is the destructor declaration

	// Functions to be defined
	bool setNewFileName(char* filePath);
	char* getParamString(const char* paraName);

	// Members
	std::ifstream paramFileStream;
	char* paramBuffer;
	char* filePathString;
};

#endif /* CPARAMREADER_H_ */
