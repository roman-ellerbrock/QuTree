//
// Created by Roman Ellerbrock on 1/21/22.
//

#include <UnitTest++/UnitTest++.h>
#include "Tree/LeafArray.h"
#include "Tree/Tree.h"
#include "Tree/SubTree.h"

SUITE (SubTree) {

	class SomeTree {
	public:
		SomeTree() {
			string file("1	-2\n"
						"	 	1	-2\n"
						"	 		1	-1\n"
						" 				1	0	0\n"
						"	 		1	-1\n"
						"	 			1	0	1\n"
						"	 	1	-2\n"
						"	 		1	-1\n"
						"	 			1	0	2\n"
						"	 		1	-1\n"
						"	 			1	0	3\n"
						"1.	0.	0.	1.\n"
						"1.	0.	0.	1.\n"
						"1.	0.	0.	1.\n"
						"1.	0.	0.	1.\n");
			stringstream ss(file);
			tree_ = Tree(ss);
		}

		~SomeTree() = default;

		Tree tree_;
	};

	TEST_FIXTURE (SomeTree, init) {
		SubTree subtree;
		subtree.initialize(tree_);

		size_t i = 0;
		for (const Node* node : subtree) {
			const Node& other = tree_.nodes()[i++];
			CHECK_EQUAL(&other, node);
		}
	}

	TEST_FIXTURE (SomeTree, initPartial) {
		vector<size_t> idx = {0, 1};
		SubTree stree(tree_, idx);
		cout << "final stree:\n";
		stree.print();
	}
}

