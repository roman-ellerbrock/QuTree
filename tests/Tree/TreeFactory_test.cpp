//
// Created by Roman Ellerbrock on 12/3/21.
//

#include <gtest/gtest.h>
#include "Tree/TreeFactory.h"

TEST(TreeFactory, binary)
{
	string file("1	-2\n"
				"	 	2	-2\n"
				"	 		2	-1\n"
				" 				3	6	0\n"
				"	 		2	-1\n"
				"	 			3	6	1\n"
				"	 	2	-2\n"
				"	 		2	-1\n"
				"	 			3	6	2\n"
				"	 		2	-1\n"
				"	 			3	6	3\n"
				"1.	0.	0.	1.\n"
				"1.	0.	0.	1.\n"
				"1.	0.	0.	1.\n"
				"1.	0.	0.	1.\n");
	stringstream is(file);
	Tree expected(is);
	Tree actual = balancedTree(4, 3, 2);
	CHECK_EQUAL(expected, actual);
}

TEST(TreeFactory, closeToBalanced)
{
	string file("1	-2\n"
				"	 	2	-2\n"
				"	 		2	-1\n"
				" 				3	6	0\n"
				"	 		2	-1\n"
				"	 			3	6	1\n"
				" 		2	-1\n"
				" 			3	6	2\n"
				"1.	0.	0.	1.\n"
				"1.	0.	0.	1.\n"
				"1.	0.	0.	1.\n");
	stringstream is(file);
	Tree expected(is);
	Tree actual = balancedTree(3, 3, 2);
	CHECK_EQUAL(expected, actual);
}
