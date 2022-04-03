//
// Created by Roman Ellerbrock on 4/24/20.
//

#ifndef STRING_EXT_H
#define STRING_EXT_H
#include <string>

namespace std {

	/// Check whether string is a natural number
	/// https://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
	bool is_number(const std::string& s);

}

#endif //STRING_EXT_H
