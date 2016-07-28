#ifndef STATE_HPP_
#define STATE_HPP_

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "cpp_utils.hpp"

namespace PLISCuda {

	class State {

		unsigned int n;    // number of sub-volumes
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
			for (uint i = 0; i < n * s; i++) {
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
			unsigned int temp = first.n;
			first.n = second.n;
			second.n = temp;

			temp = first.s;
			first.s = second.s;
			second.s = temp;

			std::swap(first.state, second.state);
		}

		// ---------- Getters ----------
		unsigned int getN() const
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
		void setN(unsigned int n)
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

	inline std::ostream& operator<<(std::ostream& os, State& s)
	{
		os << "--- State --- " << "\n";
		os << "\t" << "Subvolumes number: " << s.getN() << "\n";
		os << "\t" << "Species    number: " << s.getS() << "\n";
		os << "\t" << "State: " << "\n";
		int size = s.getN() * s.getS();
		for (int i = 0; i < size; i++) {
			os << s.getState()[i] << " ";
			if ((i + 1) % s.getN() == 0)
				os << "\n";
		}

		return os;
	}

	inline std::istream& operator>>(std::istream& is, State& state)
	{
		// parse the first line (subvolumes: s)
		std::string subvolumes;
		std::getline(is, subvolumes, ':');
		std::getline(is, subvolumes);

		int n;
		try {
			n = std::stoi(subvolumes);
		} catch (std::invalid_argument &) {
			is.setstate(std::ios::failbit);
			throw std::invalid_argument("Invalid subvolume number : " + subvolumes);
		}

		if (n < 1) {
			is.setstate(std::ios::failbit);
			throw std::invalid_argument("Invalid: n < 1");
		}

		// parse the second line (species: s)
		std::string species;
		std::getline(is, species, ':');
		std::getline(is, species);

		int s;
		try {
			s = std::stoi(species);
		} catch (std::invalid_argument &) {
			is.setstate(std::ios::failbit);
			throw std::invalid_argument("Invalid species number : " + species);
		}

		if (s < 1) {
			is.setstate(std::ios::failbit);
			throw std::invalid_argument("Invalid: s < 1");
		}

		// set r and s
		state.setN(n);
		state.setS(s);

		// allocate arrays
		state.setState(new int[n * s]());

		// eat newline
		is.get();

		// loop: parse the sub-volumes lines
		unsigned int counter = 0;
		std::string subvolume;
		std::string subvolume_line;

		while (std::getline(is, subvolume, ':')) {

			// check if subvolume number is what we expect
			try {
				counter = std::stoi(subvolume);
			} catch (std::invalid_argument &) {
				is.setstate(std::ios::failbit);
				throw std::invalid_argument("Invalid subvolume number : " + subvolume);
			}

			if (counter >= state.getN()) {
				is.setstate(std::ios::failbit);
				throw std::invalid_argument("Unexpected subvolume number. "
											"Got: " + std::to_string(counter));
			}

			is.get();    // eat space after ':'

			// get rest of the line (everything after ':') and split it
			std::getline(is, subvolume_line);

			std::vector<std::string> coeffs = split(subvolume_line, ' ');

			if (coeffs.size() != (uint)s) {
				throw std::invalid_argument("Parse of subvolume line failed "
											"(subvolume " + subvolume_line + ") has wrong length");
			}

			for (int i = 0; i < s; i++) {

				// first try to parse the number
				int c;
				try {
					c = std::stoi(coeffs.at(i));
				} catch (std::invalid_argument &) {
					is.setstate(std::ios::failbit);
					throw std::invalid_argument("parse of species count number failed");
				}

				// check if it's valid
				if (c < 0) {
					is.setstate(std::ios::failbit);
					throw std::invalid_argument("Invalid: species count < 0");
				}

				state.getState()[counter + n * i] = c;

			}

		}

		return is;
	}

}

#endif /* STATE_HPP_ */
