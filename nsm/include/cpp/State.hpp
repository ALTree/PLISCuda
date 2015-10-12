#ifndef STATE_HPP_
#define STATE_HPP_

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "cpp_utils.hpp"

namespace NSMCuda {

class State {

	int n;    // number of sub-volumes
	int s;    // number of species

	int * state;

public:

	// ---------- Constructor ----------
	State()
			: n(0), s(0), state(NULL)
	{
	}

	// ---------- Copy Constructor ----------
	State(const State& other)
			: n(other.n), s(other.s), state(new int[n * s])
	{
		for (int i = 0; i < n * s; i++) {
			state[i] = other.state[i];
		}
	}

	// ---------- Copy Assignment ----------
	State& operator=(State& other)
	{
		swap(*this, other);
		return *this;
	}

	// ---------- Destructor ----------
	~State()
	{
		delete[] state;
	}

	friend void swap(State& first, State& second)
	{
		int temp = first.n;
		first.n = second.n;
		second.n = temp;

		temp = first.s;
		first.s = second.s;
		second.s = temp;

		std::swap(first.state, second.state);
	}

	// ---------- Getters ----------
	int getN() const
	{
		return n;
	}

	int getS() const
	{
		return s;
	}

	int * getState() const
	{
		return state;
	}

	// ---------- Setters ----------
	void setN(int n)
	{
		this->n = n;
	}

	void setS(int s)
	{
		this->s = s;
	}

	void setState(int * state)
	{
		this->state = state;
	}

};

std::ostream& operator<<(std::ostream& os, State& s);
std::istream& operator>>(std::istream& is, State& state);

}

#endif /* STATE_HPP_ */
