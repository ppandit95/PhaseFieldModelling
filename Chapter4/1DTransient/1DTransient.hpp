/*
 * 1DTransient.hpp
 *
 *  Created on: 14-Apr-2021
 *      Author: pushkar
 */

#ifndef 1DTRANSIENT_HPP_
#define 1DTRANSIENT_HPP_

class 1DTransient{
private:
	unsigned int N=128;
	double dx = 1;
	double mu = 1.0;
	double dt = 0.2;
	double Ti = 1.0;
	double T[N],temp[N];
public:

};



#endif /* 1DTRANSIENT_HPP_ */
