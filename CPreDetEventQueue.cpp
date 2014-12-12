/*
 * CPreDetEventQueue.cpp
 *
 *  Created on: 11 Dec 2014
 *      Author: jtrusc
 */


#include "CPreDetEventQueue.h"

#include <cstring>

CPreDetEventQueue::CPreDetEventQueue()
{
	nEvents = arrayLength = 0;
	eventArray = NULL;
	expandArray();
}

CPreDetEventQueue::~CPreDetEventQueue()
{
	for(int i=0;i<arrayLength;i++)
		delete eventArray[i];

	delete[] eventArray;
}

// expand the current array length when more events are needed.
void CPreDetEventQueue::expandArray()
{
	// new pointer array.
	Event** temp = new Event*[arrayLength+PDEA_CHUNK];

	if(arrayLength>0)
		memcpy(temp, eventArray,arrayLength*sizeof(Event*));

	for(int i=0;i<PDEA_CHUNK;i++)
	{
		temp[arrayLength+i] = new Event;
		temp[arrayLength+i]->type = UNUSED_EVENT;
		temp[arrayLength+i]->time = 0.0; // is this the right default value??
	}

	delete[] eventArray;
	eventArray = temp;

	arrayLength+=PDEA_CHUNK;
}

// Take a copy of the next event in current and shift the next event to the end of the queue.
void CPreDetEventQueue::popEvent(Event& current)
{
	if(nEvents==0)
	{
		// no events on the queue. You shouldn't be asking, I think.
		current.type = UNUSED_EVENT;
		current.time = 0.0; // A better default value??
	}
	else
	{
		// copy the next value.
		current.type = eventArray[0]->type;
		current.time = eventArray[0]->time;
		Event* temp = eventArray[0];
		if(nEvents>1)
			memmove(eventArray, eventArray+1,(nEvents-1)*sizeof(Event*)); // conditional 'cos don't know if you can move zero bytes.

		// stick the old event to the back of the used events and NULL'ed it.
		temp->type = UNUSED_EVENT;
		eventArray[nEvents-1] = temp; // put it where the last assigned element used to be.
		nEvents--;
	}
}

void CPreDetEventQueue::addEvent(int newType, double newTime)
{
	// if the array's not long enough, expand it.
	if(nEvents==arrayLength)
		expandArray();

	// wherever it goes in the array the pointer is going to come from nEvents.
	Event* temp = eventArray[nEvents];
	temp->type = newType;
	temp->time = newTime;

	// if no events, no shuffling needed.
	if(nEvents==0)
	{
		nEvents++;
		return;
	}

	// where does the new element go...
	int insert;
	if(eventArray[0]->time>newTime)
	{
		insert = 0;  // at the start.
	} else
	{
		// bisect to find it...
		int bottom = 0;
		int top = nEvents;
		while(top-bottom>1)
		{
			int mid = (top+bottom)/2;
			if(eventArray[mid]->time<newTime)
				bottom = mid;
			else
				top = mid;
		}
		insert = top;
	}

	if(insert!=nEvents)
		memmove(eventArray+(insert+1),eventArray+insert,(nEvents-insert)*sizeof(Event*));

	eventArray[insert] = temp;

	nEvents++; // added an event.
}
