//
// Created by Roman Ellerbrock on 2/25/22.
//

#include <gtest/gtest.h>
#include "Tree/sparse_vector.h"

using namespace std;
TEST(sparse_vector, create)
{
	sparse_vector<int> vec;
	for (size_t i = 0; i < 4; ++i)
	{
		vec.push_back(2 * i, 3 * i);
	}

	for (size_t i = 0; i < vec.objects_.size(); ++i)
	{
		EXPECT_EQ(3 * i, vec[2 * i]);
	}
}
