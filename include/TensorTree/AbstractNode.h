#pragma once
#include "stdafx.h"

class AbstractNode {
public:
	AbstractNode() = default;
	virtual ~AbstractNode() = default;

	virtual size_t nTotalNodes() const = 0;
	virtual size_t nNodes() const = 0;
	virtual size_t nLeaves() const = 0;

	virtual AbstractNode *nextNode() = 0;
	virtual AbstractNode *nextNodeManthe() = 0;

	virtual void info(ostream& os) const = 0;
	virtual void Write(ostream& os) const = 0;

	virtual int NodeType() const = 0;

protected:
};

