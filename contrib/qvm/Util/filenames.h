//
// Created by Roman Ellerbrock on 2/18/21.
//

#ifndef FILENAMES_H
#define FILENAMES_H
#include "Core/stdafx.h"
#include <iostream>
#include <fstream>

bool fileExists(const string& name);

string filename(const string& name, const string& ending = "");

#endif //FILENAMES_H
