//
// Created by Roman Ellerbrock on 2/18/21.
//

#include "filenames.h"

bool fileExists(const string& name) {
	ifstream f(name.c_str());
	return f.good();
}

string filename(const string& name, const string& ending) {
	if (!fileExists(name + ending)) { return name + ending; }
	size_t i = 1;
	while (true) {
		string test = name + "." + to_string(i++) + ending;
		if (!fileExists(test)) {
			return test;
		}
	}
}