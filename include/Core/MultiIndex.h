#pragma once
#include "stdafx.h"

class MultiIndex
{
public:
	MultiIndex(vector<size_t>& bounderies_)
	{
		bounderies = bounderies_;
	  length = bounderies.size();

	  //used for fast linearization
	  linear = 0;

	  //Catch special case
	  if(length == 0) 
	  {
		  bounderies.push_back(1);
		  currentState.push_back(0);
		  length = 1;
	  }

		currentState.clear();
	  for(int i = 0; i < length; i++){
		  currentState.push_back(0);
		}

		calculateTotal();
	}

	MultiIndex(size_t size, size_t length_)
	{
		length = length_;
		
		bounderies.clear(); 
		currentState.clear();
		
	  //used for fast linearization
	  linear = 0;
	  total = 1;

	  //Catch special case
	  if(length == 0) 
	  {
		  bounderies.push_back(1);
		  currentState.push_back(0);
		  length =1;
	  }
	  else
	  {
		  for(int i = 0; i < length; i++)
		  {
			  bounderies.push_back(size);
			  total *= size;
			  currentState.push_back(0);
		  }
	  }
	}

	MultiIndex(const MultiIndex& multi)
	{
	  length = multi.length;
	  total = multi.total;

	  for(int i = 0; i < length; i++)
	  {
		  bounderies.push_back(multi.bounderies[i]);
		  currentState.push_back(0);
	  }
	
	  linear = 0;
	}

	void Zero()
	{
		for(int i = 0; i < length; i++){
			currentState[i] = 0;
		}

		linear = 0;
	}

	void Max()
	{
		for(int i = 0; i < length; i++)
			currentState[i] = bounderies[i]-1;

		linear = linearize();
	}

	void print() const {
		for(int i = 0; i < length; i++)
		{
			cout << currentState[i] << " ";
		}
		cout<< endl;
		for(int i = 0; i < length; i++){
			cout << bounderies[i] << " ";
		}
		cout<< endl;	
	}

	void operator++(int a){
		for(int i = 0; i < length; i++)
		{
			currentState[i]++;
			if(currentState[i] < bounderies[i] || i == length - 1) break;

			currentState[i] = 0;
		}

		linear++;
	}

	void operator--(int a){
		for(int i = length-1; i >= 0; i--)
		{
			currentState[i]--;
			//Thats a hack if size_t goes below zero it becomes very big
			if(currentState[i] < bounderies[i]-1 || i == 0) break;

			currentState[i] = bounderies[i]-1;
		}

		linear--;
	}

	size_t operator()(size_t a) const
	{
		assert(a < length);
		return currentState[a];
	}

	size_t operator[](size_t a) const
	{
		return currentState[a];
	}
	
	bool below() const
	{
		return currentState[length-1] < bounderies[length-1];
	}
	
	bool above() const
	{
		//Thats a hack if size_t goes below zero it becomes very big
		return currentState[0] < bounderies[0];
	}

	MultiIndex operator+(MultiIndex& multi)
	{
		vector<size_t> newVector = bounderies;
		for(int i = 0; i < multi.length; i++)
			newVector.push_back(multi.bounderies[i]);
		return MultiIndex(newVector);
	}
	
	void erase(size_t position)
	{
		bounderies.erase(bounderies.begin()+position);
		currentState.erase(currentState.begin()+position);
		length--;
		calculateTotal();
	}
	
  void erase(size_t startPosition, size_t endPosition)
	{
		bounderies.erase(bounderies.begin()+startPosition,
		                 bounderies.begin()+endPosition+1);
		currentState.erase(currentState.begin()+startPosition,
		                   currentState.begin()+endPosition+1);
		length -= endPosition - startPosition + 1;
		calculateTotal();
	}
	
  void pop_back()
	{
		bounderies.pop_back();
		currentState.pop_back();
		length--;
		calculateTotal();
	}
	
  void replace(size_t position, size_t newBoundary)
	{
		bounderies[position] = newBoundary;
		calculateTotal();
	}

	size_t linearize()
	{
		size_t out = 0;
		size_t factor = 1;
		for(size_t i = 0; i < length; i++)
		{
			out += currentState[i]*factor;
			factor *= bounderies[i];
		}
		return out;
	}
	
	size_t fastlinearize()
	{
		return linear;
	}

	size_t Before(size_t a)
	{
		size_t factor = 1;
		for(size_t i = 0; i < a; i++)
		{
			factor *= bounderies[i];
		}
		return factor;
	}

	size_t Active(size_t a)
	{
		return bounderies[a];
	}
	
	size_t After(size_t a)
	{
		size_t factor = 1;
		for(size_t i = a+1; i < length; i++)
		{
			factor *= bounderies[i];
		}
		return factor;
	}

	size_t dimtot() const
	{
		return total;
	}

protected:
	void calculateTotal()
	{
		total = 1;
		for(size_t i = 0; i < length; i++)
			total *= bounderies[i];
	}
	
	vector<size_t> currentState;
	vector<size_t> bounderies;
	size_t length;
	size_t linear;
	size_t total;
};
