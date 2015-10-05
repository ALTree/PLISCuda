#ifndef REACTIONS_HPP_
#define REACTIONS_HPP_

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "utils.hpp"

namespace NSMCuda {

class Reactions {

	int r;    // number of reactions
	int s;    // number of species

	int * reactants;
	int * products;

public:

	// ---------- Constructor ----------
	Reactions()
			: r(0), s(0), reactants(NULL), products(NULL)
	{

	}

	// ---------- Copy Constructor ----------
	Reactions(const Reactions& other)
			: r(other.r), s(other.s), reactants(new int[r * s]), products(new int[r * s])
	{
		for (int i = 0; i < r * s; i++) {
			reactants[i] = other.reactants[i];
			products[i] = other.products[i];
		}
	}

	// ---------- Copy Assignment ----------
	Reactions& operator=(Reactions& other)
	{
		swap(*this, other);
		return *this;
	}

	// ---------- Destructor ----------
	~Reactions()
	{
		delete[] reactants;
		delete[] products;
	}

	friend void swap(Reactions& first, Reactions& second)
	{
		int temp = first.r;
		first.r = second.r;
		second.r = temp;

		temp = first.s;
		first.s = second.s;
		second.s = temp;

		std::swap(first.reactants, second.reactants);
		std::swap(first.products, second.products);
	}

	// ---------- Getters ----------
	int getR() const
	{
		return r;
	}

	int getS() const
	{
		return s;
	}

	int * getReactants() const
	{
		return reactants;
	}

	int * getProducts() const
	{
		return products;
	}

	// ---------- Setters ----------
	void setR(int r)
	{
		this->r = r;
	}

	void setS(int s)
	{
		this->s = s;
	}

	void setReactants(int * reactants)
	{
		this->reactants = reactants;
	}

	void setProducts(int * products)
	{
		this->products = products;
	}

};

std::ostream& operator<<(std::ostream& os, Reactions& r)
{
	os << "--- Reactions --- " << "\n";
	os << "\t" << "Reactions number: " << r.getR() << "\n";
	os << "\t" << "Species   number: " << r.getS() << "\n";
	os << "\t" << "Reactants: " << "\n\t\t";
	int size = r.getR() * r.getS();
	for (int i = 0; i < size; i++) {
		os << r.getReactants()[i] << " ";
	}

	os << "\n\t" << "Products : " << "\n\t\t";
	for (int i = 0; i < size; i++) {
		os << r.getProducts()[i] << " ";
	}

	return os;
}

std::istream& operator>>(std::istream& is, Reactions& rs)
{
	// parse the first line (reactions: r)
	std::string reactions;
	std::getline(is, reactions, ':');
	std::getline(is, reactions);

	int r;
	try {
		r = std::stoi(reactions);
	} catch (std::invalid_argument &) {
		is.setstate(std::ios::failbit);
		throw std::invalid_argument("Invalid reactions number : " + reactions);
	}

	if (r < 1) {
		is.setstate(std::ios::failbit);
		throw std::invalid_argument("Invalid: r < 1");
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
	rs.setR(r);
	rs.setS(s);

	// allocate arrays
	rs.setReactants(new int[r * s]);
	rs.setProducts(new int[r * s]);

	// eat newline
	is.get();

	// loop: parse the reactions lines
	std::string reaction_line;
	int current_reaction = 0;

	while (std::getline(is, reaction_line)) {

		std::vector<std::string> coeffs = split(reaction_line, ' ');
		if (coeffs.size() != (2 * s + 1)) {
			throw std::invalid_argument("Parse of reactions line failed "
					"(reaction " + std::to_string(current_reaction) + ")");
		}

		// get reactans
		for (int i = 0; i < s; i++) {
			try {
				int coeff = std::stoi(coeffs.at(i));
				rs.getReactants()[current_reaction + r * i] = coeff;
			} catch (std::invalid_argument &) {
				is.setstate(std::ios::failbit);
				throw std::invalid_argument("parse of stechiometric coeff. failed "
						"(reaction " + std::to_string(current_reaction) + ")");
			}
		}

		// get products
		for (int i = 0; i < s; i++) {
			try {
				int coeff = std::stoi(coeffs.at(s + i + 1));
				rs.getProducts()[current_reaction + r * i] = coeff;
			} catch (std::invalid_argument &) {
				is.setstate(std::ios::failbit);
				throw std::invalid_argument("parse of stechiometric coeff. failed"
						"(reaction " + std::to_string(current_reaction));
			}
		}

		current_reaction++;

	}

	return is;

}

}

#endif /* REACTIONS_HPP_ */
