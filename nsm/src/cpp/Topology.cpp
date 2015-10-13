#include "../../include/cpp/Topology.hpp"

namespace NSMCuda {

std::ostream& operator<<(std::ostream& os, Topology& t)
{
	os << "--- Topology --- " << "\n";
	os << "\t" << "Subvolumes number: " << t.getN() << "\n";
	os << "\t" << "Neighbours number: " << t.getNeighboursLength() << "\n";
	os << "\t" << "Neighbours: " << "\n\t\t";
	for (int i = 0; i < t.getNeighboursLength(); i++) {
		os << t.getNeighboursArray()[i] << " ";
	}
	os << "\n";
	return os;
}

std::istream& operator>>(std::istream& is, Topology& t)
{
	// eat the first line
	std::string subvolumes;
	std::getline(is, subvolumes, ':');
	std::getline(is, subvolumes);

	// try parse-to-integer on the part after ':'
	int n;
	try {
		n = std::stoi(subvolumes);
	} catch (std::invalid_argument &) {
		is.setstate(std::ios::failbit);
		throw std::invalid_argument("Invalid subvolumes number : " + subvolumes);
	}

	// check if n > 0
	if (n < 1) {
		is.setstate(std::ios::failbit);
		throw std::invalid_argument("Invalid: n < 1");
	}

	// allocate Topology arrays
	t.setN(n);
	t.setNeighboursArray(new int[6 * n]);

	// eat newline after first line
	is.get();

	// loop: parse the "subvolume: {neighbours list}" lines
	int counter = 0;
	std::string subvolume;
	std::string neighbours_line;

	while (std::getline(is, subvolume, ':')) {

		// check if subvolume number is what we expect
		try {
			if (counter == std::stoi(subvolume)) {
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
			int sub;
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
			t.getNeighboursArray()[(counter - 1)*6 + i] = sub;
			i++;
		}


		// fill others with -1
		while (i < 6) {
			t.getNeighboursArray()[(counter - 1)*6 + i] = -1;
			i++;
		}
	}

	// no more lines (and no neighbours), we're done.

	return is;
}

}
