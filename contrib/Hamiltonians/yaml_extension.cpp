//
// Created by Roman Ellerbrock on 12/31/22.
//
#include "yaml_extension.h"

template <typename T>
T evaluate(const YAML::Node& node, const std::string& key) {
	T val;
	if (auto par = node[key]) {
		return par.as<T>();
	} else {
		std::cerr << "Did not specify key '" << key << "' in yaml node " << node << std::endl;
		exit(3);
	}
}

template size_t evaluate<size_t>(const YAML::Node& node, const std::string& key);
template double evaluate<double>(const YAML::Node& node, const std::string& key);
template char evaluate<char>(const YAML::Node& node, const std::string& key);
template bool evaluate<bool>(const YAML::Node& node, const std::string& key);
template std::string evaluate<std::string>(const YAML::Node& node, const std::string& key);

template <typename T>
T evaluate(const YAML::Node& node, const std::string& key, T default_) {
	T val;
	if (auto par = node[key]) {
		return par.as<T>();
	} else {
		return default_;
	}
}

template <typename T>
T evaluate(const YAML::Node& node, const std::string& key, T default_);

template size_t evaluate(const YAML::Node& node, const std::string& key, size_t default_);
template char evaluate(const YAML::Node& node, const std::string& key, char default_);
template bool evaluate(const YAML::Node& node, const std::string& key, bool default_);
template double evaluate(const YAML::Node& node, const std::string& key, double default_);
template std::string evaluate(const YAML::Node& node, const std::string& key, std::string default_);

/// unable to instantiate the correct templates above, so using this heler instead
template <typename T>
T instantiate_evaluate() {
	YAML::Node node;
	evaluate<std::string>(node, "string", "arg");
	evaluate<std::string>(node, "string");
	evaluate<size_t>(node, "string");
	evaluate<size_t>(node, "string", 20);
	double x = 1.;
	return x;
}

template double instantiate_evaluate<double>();
