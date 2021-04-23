/*
 * transient.hpp
 *
 *  Created on: 16-Apr-2021
 *      Author: pushkar
 */

#ifndef TRANSIENT_HPP_
#define TRANSIENT_HPP_
#include"structures.hpp"
class TransientProblem{
private:
	Params parameters;
public:
	TransientProblem();
	~TransientProblem();
	Params GetParameters() const;
	void initialize_field();
	void Output() const;
	void setup_periodic_BC();
	void Jacobi_iteration();
	void Output(unsigned int t) const;
	void retransfer();
	void solve();

};



#endif /* TRANSIENT_HPP_ */
