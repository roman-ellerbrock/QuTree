//
// Created by Roman Ellerbrock on 12/3/21.
//

#include <UnitTest++/UnitTest++.h>
#include "Tree/TreeFactory.h"


SUITE(TreeFactory) {

	TEST(binary) {
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

	TEST(closeToBalanced) {
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

}