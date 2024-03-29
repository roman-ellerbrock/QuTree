#pragma once
#include "Core/stdafx.h"

class NodePosition: public vector<size_t>
	/**
	 * \class NodePosition
	 * \ingroup TTBasis
	 * \brief This class manages the position information in a tree, mainly for I/O purposes.
	 */
{
public:
	NodePosition() {
		push_back(0);
	}

	~NodePosition() = default;

	friend NodePosition operator*(NodePosition p, size_t k);

	friend NodePosition operator*(NodePosition p, NodePosition q);

//	void push_back(int parent) { path_.push_back(parent); }

	void info(ostream& os = cout, bool print_layer = false) const;

	size_t childIdx() const;

//	size_t layer() const { return layer_; }
	size_t layer() const { return size() - 1; }

protected:
//	vector<size_t> path_;
//	size_t layer_;
};

