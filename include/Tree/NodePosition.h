#pragma once
#include "stdafx.h"

class NodePosition: public list<size_t>
	/**
	 * \class NodePosition
	 * \ingroup TTBasis
	 * \brief This class manages the position information in a tree, mainly for I/O purposes.
	 */
{
public:
	NodePosition() = default;
	~NodePosition() = default;

	explicit NodePosition(const initializer_list<size_t> p)
		: list<size_t>(p) {}

	friend NodePosition operator*(NodePosition p, size_t k);
	friend NodePosition operator*(NodePosition p, NodePosition q);
	void info(ostream& os = cout, bool print_layer = false) const;
	size_t childIdx() const;

	size_t layer() const { return size(); }
};

bool operator==(NodePosition p, NodePosition q);
ostream& operator<<(ostream& stream, const NodePosition& p);
