/*
 * File Name: CPreDetEventQueue.h
 *
 *  Created on: 14 April 2015
 *  Author: Jie Yang
 *
 */

#ifndef CPREDETEVENTQUEUE_H_
#define CPREDETEVENTQUEUE_H_

#include <cstddef> // Sometimes NULL isn't declared?!

struct Event
{
	int type;
	double time;
	void* subject; // Pointer to subject of event if needed
};

#define PreDetEvent_CHUNK 20
#define UNUSED_EVENT 0
#define RATE_EVENT 1
#define HOST_DEATH 2
#define WORM_BIRTH_DEATH 3
#define WORM_FREELIVING 4
#define TREATMENT_EVENT 5
#define TERMINATE 6
#define DEBUG_EVENT 7
#define SURVEY_EVENT 8

class CPreDetEventQueue
{
public:
CPreDetEventQueue();
	virtual ~CPreDetEventQueue();
	void popEvent(Event& current);
	void addEvent(int newType, double newTime, void* currentSubject = NULL);
private:
	void expandArray();

	// Members
	int nEvents, arrayLength;
	Event** eventArray;
};

#endif /* CPREDETEVENTQUEUE_H_ */
