#ifndef TOPOLOGY_HPP_
#define TOPOLOGY_HPP_

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "cpp_utils.hpp"

namespace NSMCuda {

class Topology {

	int n;    // number of sub-volumes

	int * neighbours_array;
	int neighbours_length;

public:

	// ---------- Constructor ----------
	Topology()
			: n(0), neighbours_array(NULL), neighbours_length(0)
	{

	}

	// ---------- Copy Constructor ----------
	Topology(const Topology& other)
			: n(other.n), neighbours_array(new int[n * 6]), neighbours_length(6 * n)
	{
		for (int i = 0; i < 6 * n; i++) {
			neighbours_array[i] = other.neighbours_array[i];
		}
	}

	// ---------- Copy Assignment ----------
	Topology& operator=(Topology& other)
	{
		swap(*this, other);
		return *this;
	}

	// ---------- Destructor ----------
	~Topology()
	{
		delete[] neighbours_array;
	}

	friend void swap(Topology& first, Topology& second)
	{
		int temp = first.n;
		first.n = second.n;
		second.n = temp;

		std::swap(first.neighbours_array, second.neighbours_array);

		temp = first.neighbours_length;
		first.neighbours_length = second.neighbours_length;
		second.neighbours_length = temp;
	}

	// ---------- Getters ----------
	int getN() const
	{
		return n;
	}

	int * getNeighboursArray() const
	{
		return neighbours_array;
	}

	int getNeighboursLength() const
	{
		return neighbours_length;
	}

	// ---------- Setters ----------
	void setN(int n)
	{
		this->n = n;
		this->neighbours_length = 6*n;
	}

	void setNeighboursArray(int * neighboursArray)
	{
		neighbours_array = neighboursArray;
	}

};

std::ostream& operator<<(std::ostream& os, NSMCuda::Topology& t);
std::istream& operator>>(std::istream& is, Topology& t);

}

#endif /* TOPOLOGY_HPP_ */
