#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <string>
#include <sstream>
#include <vector>

namespace NSMCuda {

// http://stackoverflow.com/a/236803/1146303
// inline to avoid 'multiple definition' error:
//		http://stackoverflow.com/a/14425299/1146303
// TODO: remove
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> split(const std::string &s, char delim);

}

#endif /* UTILS_HPP_ */
