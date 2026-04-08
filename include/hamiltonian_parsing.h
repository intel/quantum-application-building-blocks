//===--- hamiltonian_parsing.h ----------------------------------*- C++ -*-===//
//
//===----------------------------------------------------------------------===//
//
// Copyright (C) 2026 Intel Corporation.
//
// This software and the related documents are Intel copyrighted materials, and
// your use of them is governed by the express license under which they were
// provided to you ("License"). Unless the License provides otherwise, you may
// not use, modify, copy, publish, distribute, disclose or transmit this
// software or the related documents without Intel's prior written permission.
//
// This software and the related documents are provided as is, with no express
// or implied warranties, other than those that are expressly stated in the
// License.
//===----------------------------------------------------------------------===//
///
/// \file hamiltonian_parsing.h
/// \brief A library of functions and classes to represent Hamiltonians as
///  Pauli strings.
///
/// A set of tools to represent Hamiltonians as Pauli strings and read formatted
/// files containing one or more terms of Pauli strings.
///
//===----------------------------------------------------------------------===//

#ifndef HAM_PARSE_H
#define HAM_PARSE_H

// C++ standard library
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// @brief A class to represent one term out of a linear combination of Pauli
///        strings that make up a physics Hamiltonian
struct HamiltonianTerm {
    double coeff;  /*!< numerical coefficient or weight of this term within a
                        larger expansion */
    std::string pauli; /*!< string representing a set of Pauli gates applied to
                            a set of qubits of equal length */
};

/// @brief Load Hamiltonian from file
/// @param filename The path and name of the file to read.
std::vector<HamiltonianTerm> load_hamiltonian(const std::string& filename) {

    std::vector<HamiltonianTerm> terms;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Check if file exists and is readable." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    unsigned line_number = 0;

    while (std::getline(file, line)) {
        line_number++;
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

        if (line.empty()) continue;

        size_t star_pos = line.find('*');
        if (star_pos == std::string::npos) {
            std::cerr << "Warning: Invalid format at line " << line_number 
                      << " (missing '*'): " << line << std::endl;
            continue;
        }

        try {
            // coefficient
            std::string coeff_str = line.substr(0, star_pos);
            if (coeff_str[0] == '+') {
                coeff_str = coeff_str.substr(1); // Remove leading +
            }

            // Pauli string
            std::string pauli_str = line.substr(star_pos + 1);

            // validate Pauli string
            for (char c : pauli_str) {
                if (c != 'I' && c != 'X' && c != 'Y' && c != 'Z') {
                    std::cerr << "Warning: Invalid Pauli operator '" << c 
                              << "' at line " << line_number << std::endl;
                }
            }

            HamiltonianTerm term;
            term.coeff = std::stod(coeff_str);
            term.pauli = pauli_str;
            terms.push_back(term);

        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to parse line " << line_number 
                      << ": " << line << " (" << e.what() << ")" << std::endl;
        }
    }

    file.close();

    // ensure we loaded at least one term
    if (terms.empty()) {
        std::cerr << "Error: No valid Hamiltonian terms found in file: " << filename << std::endl;
        std::cerr << "Expected format: coefficient * PauliString" << std::endl;
        std::cerr << "Example: -1.0523732 * II" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return terms;
}

#endif // HAM_PARSE_H
