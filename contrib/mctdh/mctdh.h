//
// Created by Roman Ellerbrock on 2/27/20.
//

#ifndef YAML_PARSER_H
#define YAML_PARSER_H
#include "mctdh_state.h"
#include "mctdh.h"
#include "yaml-cpp/yaml.h"

namespace parser {

	mctdh_state read(const string& yaml_filename);

	mctdh_state run(const string& yaml_filename);

	Tree create_tree(const YAML::Node& node);
}

#endif //YAML_PARSER_H
