//
// Created by Roman Ellerbrock on 2/25/22.
//

#include <UnitTest++/UnitTest++.h>
#include "Tree/sparse_vector.h"


SUITE(sparse_vector) {
	using namespace std;
	TEST(create) {
		sparse_vector<int> vec;
		for (size_t i = 0; i < 4; ++i) {
			vec.push_back(2*i, 3*i);
		}

		for (size_t i = 0; i < vec.objects_.size(); ++i) {
			CHECK_EQUAL(3 * i, vec[2*i]);
		}

	}
}
