//
// Created by Roman Ellerbrock on 12/31/22.
//

#ifndef YAML_EXTENSION_H
#define YAML_EXTENSION_H
#include "yaml-cpp/yaml.h"
#include <string>
#include <iostream>

template <typename T>
T evaluate(const YAML::Node& node, const std::string& key);

template <typename T>
T evaluate(const YAML::Node& node, const std::string& key, T default_);

#endif //YAML_EXTENSION_H
