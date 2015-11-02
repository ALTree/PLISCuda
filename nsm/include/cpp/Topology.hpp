#ifndef TOPOLOGY_HPP_
#define TOPOLOGY_HPP_

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "cpp_utils.hpp"

namespace NSMCuda {

class Topology {

	unsigned int n;    // number of sub-volumes

	int * neighbours_array;

public:

	// ---------- Constructor ----------
	Topology()
			: n(0), neighbours_array(NULL)
	{

	}

	// ---------- Copy Constructor ----------
	Topology(const Topology& other)
			: n(other.n), neighbours_array(new int[n * 6])
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
		unsigned int temp = first.n;
		first.n = second.n;
		second.n = temp;

		std::swap(first.neighbours_array, second.neighbours_array);
	}

	// ---------- Getters ----------
	unsigned int getN() const
	{
		return n;
	}

	int * getNeighboursArray() const
	{
		return neighbours_array;
	}

	int getNeighboursLength() const
	{
		return 6 * n;
	}

	// ---------- Setters ----------
	void setN(int n)
	{
		this->n = n;
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
