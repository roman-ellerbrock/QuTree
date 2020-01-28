#pragma once
#include "Leaf.h"

class LinearizedLeaves
{
public:
	LinearizedLeaves() = default;
	~LinearizedLeaves() = default;

	size_t size()const { return coordinates_.size(); }

	void push_back(Leaf& phys);
	void resizeaddress(int n);
	void clear();
	const Leaf& operator[](int mode)const;
	Leaf& operator[](int mode);
	reference_wrapper<Leaf>& operator()(int mode);

protected:
	vector<reference_wrapper<Leaf>> coordinates_;
	vector<int> address_;
};

