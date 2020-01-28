#pragma once
#include "Core/stdafx.h"

class NodePosition
{
 public:
	NodePosition():layer_(0) {}
	~NodePosition() = default;

	friend NodePosition operator*(NodePosition p, int k);

	friend NodePosition operator*(NodePosition p, NodePosition q);

	void push_back(int parent) { path_.push_back(parent); }

	void info(ostream& os=cout)const;

	int ChildIdx()const;

	size_t Layer()const { return layer_; }

 protected:
	vector<int> path_;
	int layer_;
};

