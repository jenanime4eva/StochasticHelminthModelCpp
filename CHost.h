/*
 * CHost.h
 *
 *  Created on: 1 Dec 2014
 *      Author: Jie yang
 */

#ifndef CHOST_H_
#define CHOST_H_

using namespace std;

class CHost {
public:
	CHost();
	virtual ~CHost();

	// Members.
	double birthDate, deathDate;
	int femaleWorms, totalWorms;
	int freeliving;
};

#endif /* CHOST_H_ */
