/*
 * CPreDetEventQueue.h
 *
 *  Created on: 11 Dec 2014
 *      Author: jtrusc
 */

#ifndef CPREDETEVENTQUEUE_H_
#define CPREDETEVENTQUEUE_H_

#include <cstddef>  // Sometimes NULL isn't declared??!?

struct Event
{
	int type;
	double time;
	void* subject; // pointer to subject of event, if needed.
};

#define PDEA_CHUNK 20
#define UNUSED_EVENT 0
#define HOST_DEATH 1
#define TERMINATE 2
#define DEBUG_EVENT 3
#define SURVEY_EVENT 4

class CPreDetEventQueue
{
public:
	CPreDetEventQueue();
	virtual ~CPreDetEventQueue();
	void popEvent(Event& current);
	void addEvent(int newType, double newTime, void* currentSubject = NULL);

private:
	void expandArray();

	// members...
	int nEvents, arrayLength;
	Event** eventArray;
};

#endif /* CPREDETEVENTQUEUE_H_ */