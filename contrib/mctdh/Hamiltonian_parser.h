//
// Created by Roman Ellerbrock on 12/31/22.
//

#ifndef HAMILTONIAN_PARSER_H
#define HAMILTONIAN_PARSER_H
#include "Hamiltonians.h"
#include "yaml-cpp/yaml.h"

shared_ptr<Hamiltonian> read_hamiltonian(const YAML::Node& node,
	const Tree& tree);

PotentialOperator set_potential(const YAML::Node& node, const Tree& tree);

#endif //HAMILTONIAN_PARSER_H
