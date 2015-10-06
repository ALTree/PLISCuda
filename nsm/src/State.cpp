#include "../include/State.hpp"

namespace NSMCuda {

std::ostream& operator<<(std::ostream& os, State& s)
{
	os << "--- State --- " << "\n";
	os << "\t" << "Subvolumes number: " << s.getN() << "\n";
	os << "\t" << "Species    number: " << s.getS() << "\n";
	os << "\t" << "State: " << "\n\t\t";
	int size = s.getN() * s.getS();
	for (int i = 0; i < size; i++) {
		os << s.getState()[i] << " ";
		if (i % s.getS() == 0)
			os << "\n";
	}

	return os;
}

std::istream& operator>>(std::istream& is, State& state)
{
	// parse the first line (reactions: r)
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
	state.setState(new int[n * s]);

	// eat newline
	is.get();

	// loop: parse the sub-volumes lines
	int counter = 0;
	std::string subvolume;
	std::string subvolume_line;
	int current_subvolume = 0;

	while (std::getline(is, subvolume, ':')) {

		// check if subvolume number is what we expect
		try {
			if (counter == std::stoi(subvolume)) {
				// yes. write current subvolume in offset array and increment counter
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
		std::getline(is, subvolume_line);

		std::vector<std::string> coeffs = split(subvolume_line, ' ');

		if (coeffs.size() != s) {
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

			state.getState()[current_subvolume + n * i] = c;

		}

		current_subvolume++;
	}

	return is;
}

}
