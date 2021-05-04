/*
 * structures.hpp
 *
 *  Created on: 03-May-2021
 *      Author: pushkar
 */

#ifndef SRC_STRUCTURES_HPP_
#define SRC_STRUCTURES_HPP_

struct SimParams{
	unsigned int Nx;
	unsigned int Ny;
	double dx;
	double dy;
	double** con,**mu,**dfdc,**laplace_con,**laplace_dfdc,**con_temp;
	char path[256];
};
struct TimeParams{
	unsigned int nstep;
	double dtime;
	double ttime;
};
struct MatParams{
	double c0;
	double mobility;
	double grad_coef;
	double A;
	double energy;
};



#endif /* SRC_STRUCTURES_HPP_ */
