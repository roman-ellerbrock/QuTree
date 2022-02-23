#pragma once
#include "Leaf.h"
#include "Node.h"

class LeafArray
	/**
	 * \class LinearizedLeaves
	 * \ingroup Tree
	 * \brief The class holds a vector of references to the leaves in the Tree.
	 */
{
public:
	LeafArray() = default;
	~LeafArray() = default;
	explicit LeafArray(Node& root);

	[[nodiscard]] size_t size()const { return leaves_.size(); }

	void readPars(istream& is) {
		for (Leaf& leaf : leaves_) {
			leaf.readPar(is);
		}
	}

	void writePars(ostream& os) const {
		for (const Leaf& leaf : leaves_) {
			leaf.writePar(os);
		}
	}

	void push_back(Leaf& leaf);

	const Leaf& operator[](size_t mode)const;
	Leaf& operator[](size_t mode);

protected:
	vector<reference_wrapper<Leaf>> leaves_;
	vector<int> address_;
};

