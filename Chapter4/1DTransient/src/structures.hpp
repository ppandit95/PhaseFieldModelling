/*
 * structures.hpp
 *
 *  Created on: 16-Apr-2021
 *      Author: pushkar
 */

#ifndef STRUCTURES_HPP_
#define STRUCTURES_HPP_


struct Params{
	unsigned int N;
	double dx;
	double mu;
	double dt;
	double Ti;
	double *T,*temp;
};


#endif /* STRUCTURES_HPP_ */
