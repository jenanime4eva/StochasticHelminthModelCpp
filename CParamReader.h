/*
 * CParamReader.h
 *
 *  Created on: 12 Nov 2014
 *      Author: jtrusc
 *	Modified: 26 Nov 2014
 *		Author: Jie Yang
 *      Each informative line of the param file should have the form:
 *      paramNumber (tab) Value (tab) ## Comments
 *      Any line starting with '#' is completely ignored.
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

	bool setNewFileName(char* filePath);
	char* getParamString(const char* paraName);


	// Members...
	std::ifstream paramFileStream;
	char* paramBuffer;
	char* filePathString;
};

#endif /* CPARAMREADER_H_ */
