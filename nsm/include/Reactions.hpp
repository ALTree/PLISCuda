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

std::ostream& operator<<(std::ostream& os, Reactions& r);
std::istream& operator>>(std::istream& is, Reactions& rs);

}

#endif /* REACTIONS_HPP_ */
