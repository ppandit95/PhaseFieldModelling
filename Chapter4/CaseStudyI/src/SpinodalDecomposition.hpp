/*
 * SpinodalDecomposition.hpp
 *
 *  Created on: 03-May-2021
 *      Author: pushkar
 */

#ifndef SRC_SPINODALDECOMPOSITION_HPP_
#define SRC_SPINODALDECOMPOSITION_HPP_
#include"structures.hpp"

class SpinodalDecomposition{
private:
	SimParams p;
	TimeParams t;
	MatParams m;
public:
	SpinodalDecomposition();
	~SpinodalDecomposition();
	SimParams GetSimParams();
	TimeParams GetTimeParams();
	MatParams GetMatParams();
	void Initial_Profile();
	void write_output(std::string filename);
	void CalcBulkEnergy(unsigned int i,unsigned int j);
	void CalcChemicalPotential(unsigned int i,unsigned int j);
	void Set_Periodic_BCs(unsigned int i,unsigned int j,unsigned int count);
};




#endif /* SRC_SPINODALDECOMPOSITION_HPP_ */
