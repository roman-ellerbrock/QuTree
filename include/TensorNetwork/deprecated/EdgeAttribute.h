//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef EDGEATTRIBUTE_H
#define EDGEATTRIBUTE_H


template<class A>
class EdgeAttribute : public vector<A>{
public:
	A& operator[](const Edge& e) {
		const Node& x = e.Down();
		size_t address = x.Address();
		assert(address < attributes_.size());
		return attributes_[address];
	}

	const A& operator[](const Edge& e) const {
		const Node& x = e.Down();
		size_t address = x.Address();
		assert(address < attributes_.size());
		return attributes_[address];
	}

};


#endif //EDGEATTRIBUTE_H
