#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <string>
#include <sstream>
#include <vector>

namespace NSMCuda {

// http://stackoverflow.com/a/236803/1146303
// inline to avoid 'multiple definition' error:
//		http://stackoverflow.com/a/14425299/1146303
inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

inline std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

}

#endif /* UTILS_HPP_ */
