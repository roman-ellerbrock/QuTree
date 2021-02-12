#pragma once
#include "Core/stdafx.h"
#include "NodePosition.h"

/**
 * \defgroup TTBasis
 * \brief This group contains classes to construct and handle a tensor tree basis.
 *
 */

class AbstractNode {
	/**
	 * \class AbstractNode
	 * \ingroup TTBasis
	 * \brief This is an abstract node in the tree of a TTBasis.
	 */
public:
	AbstractNode() = default;
	virtual ~AbstractNode() = default;

	virtual size_t nTotalNodes() const = 0;
	virtual size_t nNodes() const = 0;
	virtual size_t nLeaves() const = 0;

	virtual AbstractNode *nextNode() = 0;
	virtual AbstractNode *nextNodeManthe() = 0;

	virtual void info(ostream& os) const = 0;
	virtual void write(ostream& os) const = 0;

	virtual int type() const = 0;
	virtual void setParent(AbstractNode* Up) = 0;
	virtual void update(const NodePosition& p) = 0;

protected:
};

