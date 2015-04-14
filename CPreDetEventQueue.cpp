/*
 * File Name: CPreDetEventQueue.cpp
 *
 *  Created on: 14 April 2015
 *  Author: Jie Yang
 *
 */

#include "CPreDetEventQueue.h"
#include <cstring>

// Class constructor
CPreDetEventQueue::CPreDetEventQueue()
{
	nEvents = arrayLength = 0;
	eventArray = NULL;
	expandArray();
}

// Class destructor
CPreDetEventQueue::~CPreDetEventQueue()
{
	for(int i=0;i<arrayLength;i++)
		delete eventArray[i];

	delete[] eventArray;
}

// Expand the current array length when more events are needed
void CPreDetEventQueue::expandArray()
{
	// New pointer array
	Event** temp = new Event*[arrayLength+PreDetEvent_CHUNK];

	if(arrayLength>0)
		memcpy(temp, eventArray,arrayLength*sizeof(Event*));

	for(int i=0;i<PreDetEvent_CHUNK;i++)
	{
		temp[arrayLength+i] = new Event;
		temp[arrayLength+i]->type = UNUSED_EVENT;
		temp[arrayLength+i]->subject = NULL;
		temp[arrayLength+i]->time = 0.0; // Is this the right default value?
	}

	delete[] eventArray;
	eventArray = temp;

	arrayLength += PreDetEvent_CHUNK;
}

// Take a copy of the next event in current and shift the next event to the end of the queue
void CPreDetEventQueue::popEvent(Event& current)
{
	if(nEvents==0)
	{
		// No events on the queue
		current.type = UNUSED_EVENT;
		current.time = 0.0; // A better default value?
		current.subject = NULL;
	}
	else
	{
	// Copy the next value
		current.type = eventArray[0]->type;
		current.time = eventArray[0]->time;
		current.subject = eventArray[0]->subject;

		Event* temp = eventArray[0];
		if(nEvents>1)
			memmove(eventArray, eventArray+1,(nEvents-1)*sizeof(Event*)); // Conditional (don't know if you can move zero bytes)

		// Stick the old event to the back of the used events and NULL it
		temp->type = UNUSED_EVENT;
		eventArray[nEvents-1] = temp; // Put it where the last assigned element used to be
		nEvents--;
	}
}
void CPreDetEventQueue::addEvent(int newType, double newTime, void* newSubject)
{
	// If the array's not long enough expand it
	if(nEvents==arrayLength)
		expandArray();

	// Wherever it goes in the array the pointer is going to come from nEvents
	Event* temp = eventArray[nEvents];
	temp->type = newType;
	temp->time = newTime;
	temp->subject = newSubject;

	// If no events no shuffling needed
	if(nEvents==0)
	{
		nEvents++;
		return;
	}

	// Where does the new element go?
	int insert;
	if(eventArray[0]->time>newTime)
	{
		insert = 0; // At the start
	} else
	{
		// Bisect to find it
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

	nEvents++; // Added an event
}
