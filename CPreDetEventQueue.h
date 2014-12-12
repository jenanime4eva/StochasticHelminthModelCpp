/*
 * CPreDetEventQueue.h
 *
 *  Created on: 11 Dec 2014
 *      Author: jtrusc
 */

#ifndef CPREDETEVENTQUEUE_H_
#define CPREDETEVENTQUEUE_H_

struct Event
{
	int type;
	double time;
};

#define PDEA_CHUNK 5
#define UNUSED_EVENT 0
#define HOST_DEATH 1

class CPreDetEventQueue
{
public:
	CPreDetEventQueue();
	virtual ~CPreDetEventQueue();
	void popEvent(Event& current);
	void addEvent(int newType, double newTime);

private:
	void expandArray();

	// members...
	int nEvents, arrayLength;
	Event** eventArray;
};

#endif /* CPREDETEVENTQUEUE_H_ */
