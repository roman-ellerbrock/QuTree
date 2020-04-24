//
// Created by Roman Ellerbrock on 4/24/20.
//
#include "Util/string_ext.h"

namespace std {

	/// Check whether string is a natural number
	/// https://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
	bool is_number(const std::string& s)
	{
		return !s.empty() && std::find_if(s.begin(),
			s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
	}
}
