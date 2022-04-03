#pragma once
#include "stdafx.h"

class MultiIndex
{
public:
	MultiIndex(vector<size_t>& boundaries)
	{
		boundaries_ = boundaries;
	  length_ = boundaries_.size();

	  //used for fast linear_ization
	  linear_ = 0;

	  //Catch special case
	  if(length_ == 0)
	  {
		  boundaries_.push_back(1);
		  currentState_.push_back(0);
		  length_ = 1;
	  }

		currentState_.clear();
	  for(int i = 0; i < length_; i++){
		  currentState_.push_back(0);
		}

		calculateTotal();
	}

	MultiIndex(size_t size, size_t length)
	{
		length_ = length;
		
		boundaries_.clear();
		currentState_.clear();
		
	  //used for fast linear_ization
	  linear_ = 0;
	  total_ = 1;

	  //Catch special case
	  if(length_ == 0)
	  {
		  boundaries_.push_back(1);
		  currentState_.push_back(0);
		  length_ =1;
	  }
	  else
	  {
		  for(int i = 0; i < length_; i++)
		  {
			  boundaries_.push_back(size);
			  total_ *= size;
			  currentState_.push_back(0);
		  }
	  }
	}

	MultiIndex(const MultiIndex& multi)
	{
	  length_ = multi.length_;
	  total_ = multi.total_;

	  for(int i = 0; i < length_; i++)
	  {
		  boundaries_.push_back(multi.boundaries_[i]);
		  currentState_.push_back(0);
	  }
	
	  linear_ = 0;
	}

	void Zero()
	{
		for(int i = 0; i < length_; i++){
			currentState_[i] = 0;
		}

		linear_ = 0;
	}

	void Max()
	{
		for(int i = 0; i < length_; i++)
			currentState_[i] = boundaries_[i]-1;

		linear_ = Linearize();
	}

	void print() const {
		for(int i = 0; i < length_; i++)
		{
			cout << currentState_[i] << " ";
		}
		cout<< endl;
		for(int i = 0; i < length_; i++){
			cout << boundaries_[i] << " ";
		}
		cout<< endl;	
	}

	void operator++(int a){
		for(int i = 0; i < length_; i++)
		{
			currentState_[i]++;
			if(currentState_[i] < boundaries_[i] || i == length_ - 1) break;

			currentState_[i] = 0;
		}

		linear_++;
	}

	void operator--(int a){
		for(int i = length_-1; i >= 0; i--)
		{
			currentState_[i]--;
			//Thats a hack if size_t goes Below zero it becomes very big
			if(currentState_[i] < boundaries_[i]-1 || i == 0) break;

			currentState_[i] = boundaries_[i]-1;
		}

		linear_--;
	}

	size_t operator()(size_t a) const
	{
		assert(a < length_);
		return currentState_[a];
	}

	size_t operator[](size_t a) const
	{
		return currentState_[a];
	}
	
	bool Below() const
	{
		return currentState_[length_-1] < boundaries_[length_-1];
	}
	
	bool After() const
	{
		//Thats a hack if size_t goes Below zero it becomes very big
		return currentState_[0] < boundaries_[0];
	}

	MultiIndex operator+(MultiIndex& multi)
	{
		vector<size_t> newVector = boundaries_;
		for(int i = 0; i < multi.length_; i++)
			newVector.push_back(multi.boundaries_[i]);
		return MultiIndex(newVector);
	}
	
	void Erase(size_t position)
	{
		boundaries_.Erase(boundaries_.begin()+position);
		currentState_.Erase(currentState_.begin()+position);
		length_--;
		calculateTotal();
	}
	
  void Erase(size_t startPosition, size_t endPosition)
	{
		boundaries_.Erase(boundaries_.begin()+startPosition,
		                 boundaries_.begin()+endPosition+1);
		currentState_.Erase(currentState_.begin()+startPosition,
		                   currentState_.begin()+endPosition+1);
		length_ -= endPosition - startPosition + 1;
		calculateTotal();
	}
	
  void pop_back()
	{
		boundaries_.pop_back();
		currentState_.pop_back();
		length_--;
		calculateTotal();
	}
	
  void Replace(size_t position, size_t newBoundary)
	{
		boundaries_[position] = newBoundary;
		calculateTotal();
	}

	size_t Linearize()
	{
		size_t out = 0;
		size_t factor = 1;
		for(size_t i = 0; i < length_; i++)
		{
			out += currentState_[i]*factor;
			factor *= boundaries_[i];
		}
		return out;
	}
	
	size_t FastLinearize()
	{
		return linear_;
	}

	size_t Before(size_t a)
	{
		size_t factor = 1;
		for(size_t i = 0; i < a; i++)
		{
			factor *= boundaries_[i];
		}
		return factor;
	}

	size_t Active(size_t a)
	{
		return boundaries_[a];
	}
	
	size_t After(size_t a)
	{
		size_t factor = 1;
		for(size_t i = a+1; i < length_; i++)
		{
			factor *= boundaries_[i];
		}
		return factor;
	}

	size_t DimTot() const
	{
		return total_;
	}

protected:
	void calculateTotal()
	{
		total_ = 1;
		for(size_t i = 0; i < length_; i++)
			total_ *= boundaries_[i];
	}
	
	vector<size_t> currentState_;
	vector<size_t> boundaries_;
	size_t length_;
	size_t linear_;
	size_t total_;
};
