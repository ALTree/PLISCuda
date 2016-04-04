#ifndef TOPOLOGY_HPP_
#define TOPOLOGY_HPP_

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "cpp_utils.hpp"

namespace PLISCuda {

class Topology {

	unsigned int n;    // number of sub-volumes
	unsigned int * neighbours_array;

public:

	// ---------- Constructor ----------
	Topology()
			: n(0), neighbours_array(NULL)
	{

	}

	// ---------- Copy Constructor ----------
	Topology(const Topology& other)
			: n(other.n), neighbours_array(new unsigned int[n * 6])
	{
		for (unsigned int i = 0; i < 6 * n; i++) {
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

	unsigned int * getNeighboursArray() const
	{
		return neighbours_array;
	}

	unsigned int getNeighboursLength() const
	{
		return 6 * n;
	}

	// ---------- Setters ----------
	void setN(unsigned int n)
	{
		this->n = n;
	}

	void setNeighboursArray(unsigned int * neighboursArray)
	{
		neighbours_array = neighboursArray;
	}

};

inline std::ostream& operator<<(std::ostream& os, Topology& t)
{
	os << "--- Topology --- " << "\n";
	os << "\t" << "Subvolumes number: " << t.getN() << "\n";
	os << "\t" << "Neighbours: " << "\n\t\t";
	for (unsigned int i = 0; i < t.getNeighboursLength(); i++) {
		os << t.getNeighboursArray()[i] << " ";
	}
	os << "\n";
	return os;
}

inline std::istream& operator>>(std::istream& is, Topology& t)
{
	// eat the first line
	std::string subvolumes;
	std::getline(is, subvolumes, ':');
	std::getline(is, subvolumes);

	// try parse-to-integer on the part after ':'
	unsigned int n;
	try {
		n = std::stoi(subvolumes);
	} catch (std::invalid_argument &) {
		is.setstate(std::ios::failbit);
		throw std::invalid_argument("Invalid subvolumes number: " + subvolumes);
	}

	// check if n > 0
	if (n < 1) {
		is.setstate(std::ios::failbit);
		throw std::invalid_argument("Invalid: n < 1");
	}

	// allocate Topology arrays
	t.setN(n);
	t.setNeighboursArray(new unsigned int[6 * n]);

	// eat newline after first line
	is.get();

	// loop: parse the "subvolume: {neighbours list}" lines
	unsigned int counter = 0;
	std::string subvolume;
	std::string neighbours_line;

	while (std::getline(is, subvolume, ':')) {

		// check if subvolume number is what we expect
		try {
			if (counter == (uint)std::stoi(subvolume)) {
				counter++;
			} else {
				// no. abort parsing.
				is.setstate(std::ios::failbit);
				throw std::invalid_argument("Unexpected subvolume number. "
						"Want: " + std::to_string(counter) + "Got: " + subvolume);
			}
		} catch (std::invalid_argument &) {
			is.setstate(std::ios::failbit);
			throw std::invalid_argument("Invalid subvolume number : " + subvolume);
		}

		is.get();    // eat space after ':'

		// get rest of the line (everything after ':') and split it
		std::getline(is, neighbours_line);
		std::vector<std::string> neighbours = split(neighbours_line, ' ');

		// loop over line elements, convert to int
		// and add to neighbours array
		int i = 0;
		for (auto& subv : neighbours) {

			// first try to parse the number
			uint sub;
			try {
				sub = std::stoi(subv);
			} catch (std::invalid_argument &) {
				is.setstate(std::ios::failbit);
				throw std::invalid_argument("parse of neighbour number failed");
			}

			// check if it's valid
			if (sub >= n) {
				is.setstate(std::ios::failbit);
				throw std::invalid_argument("Invalid: sub >= n");
			}

			// add to neighbours array and increment counter
			t.getNeighboursArray()[(counter - 1) * 6 + i] = sub;
			i++;
		}

		// fill others with -1
		unsigned int current_sbi = std::stoi(subvolume);
		while (i < 6) {
			t.getNeighboursArray()[(counter - 1) * 6 + i] = current_sbi;
			i++;
		}
	}

	// no more lines (and no neighbours), we're done.

	return is;
}

}

#endif /* TOPOLOGY_HPP_ */
