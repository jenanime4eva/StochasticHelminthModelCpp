/*
 * CHost.cpp
 *
 *  Created on: 27 Nov 2014
 *      Author: Jie Yang
 */

#include "CHost.h"
#include "CParamReader.h"
#include <iostream>
#include <string.h>

// Class Constructor
CHost::CHost()
{
	hostPopulation = NULL;
}

// Class Destructor
CHost::~CHost()
{
	if(hostPopulation!=NULL)
		delete[] hostPopulation;
}

