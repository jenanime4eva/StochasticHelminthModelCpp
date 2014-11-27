/*
 * CSimulator.cpp
 *
 *  Created on: 26 Nov 2014
 *      Author: Jie Yang
 */

#include "CSimulator.h"
#include "CParamReader.h"
#include <iostream>
#include <string.h>

// Class Constructor
CSimulator::CSimulator()
{
	results = NULL;
}

// Class Destructor
CSimulator::~CSimulator()
{
	// Delete results allocations
	if (results!=NULL)
	{
		// Delete memory
		for (int i=0 < nRepetitions;i++)
		{
			delete[] results[i];
		}

		delete[] results;
	}
}

