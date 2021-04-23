/*
 * Tests.cpp
 *
 *  Created on: 23-Apr-2021
 *      Author: pushkar
 */

#define BOOST_TEST_MODULE sample
#include "boost/test/included/unit_test.hpp"
#include "../src/transient.hpp"
#include"../src/structures.hpp"
TransientProblem	dummy;
Params  p = dummy.GetParameters();

BOOST_AUTO_TEST_CASE(init_field){
	dummy.initialize_field();
	for(int i=0;i<20;i++){
		BOOST_CHECK(p.T[(p.N/2) + i] == 1.0);
		BOOST_CHECK(p.T[(p.N/2) - i] == 1.0);
	}
}



