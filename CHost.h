/*
 * CHost.h
 *
 *  Created on: 1 Dec 2014
 *      Author: Jie yang
 */

#ifndef CHOST_H_
#define CHOST_H_

class CHost {
public:
	CHost();
	virtual ~CHost();

	// Members.
	double birthDate, deathDate;
	double femaleWorms, totalWorms;
	double wormTotalDeathRate, hostInfectionRate;
};

#endif /* CHOST_H_ */
