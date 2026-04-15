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
/// \file vqe_qk.cpp
/// \brief Demonstrate the hybrid quantum-classical algorithm
///        Variational Quantum Eigensolver using ``quantum_kernel`` functions.
///
/// Variational Quantum Eigensolver
/// Peruzzo, A., McClean, J., Shadbolt, P. et al. A variational eigenvalue
///   solver on a photonic quantum processor. Nat Commun 5, 4213 (2014).
///   https://doi.org/10.1038/ncomms5213
/// (VQE) is a hybrid quantum-classical algorithm that computes the value of a
/// cost function on the quantum hardware and optimizes the input values to find
/// an extremum of the cost function. This extremum is partly determined by the
/// form of the trial eigenfunction.
///
/// This example of VQE solves for the energy of a molecule by using the
/// Hamiltonian (formal energy equation in physics) as the cost function, and
/// a generic circuit as the trial eigenfunction, or trial wave function in
/// physics and chemistry language. The value of the cost function is optimized
/// by the SPSA method to find the rotation angles that produce the lowest
/// total energy.
///
///  Hamiltonian data obtained from HamLib library:
///  https://portal.nersc.gov/cfs/m888/dcamps/hamlib/
///    HamLib: A library of Hamiltonians for benchmarking quantum algorithms and
///      hardware
///    Nicolas P. D. Sawaya, Daniel Marti-Dafcik, Yang Ho, Daniel P Tabor, 
///      David Bernal, Alicia B Magann, Shavindra Premaratne, Pradeep Dubey,
///      Anne Matsuura, Nathan Bishop, Wibe A. de Jong, Simon Benjamin,
///      Ojas D Parekh, Norm M. Tubman, Katherine Klymko, Daan Camps
///    December 2024.
///
///  For a discussion of H2 and LiH in the spirit of this example, see
///    McArdle, S., Endo, S., Aspuru-Guzik, A., Benjamin, S.C., Yuan, X..
///      Quantum computational chemistry. Reviews of Modern Physics
///      2020;92(1):15003–15003. doi:10.1103/revmodphys.92.015003.
///
/// In this example, the circuit representing the wave function has been written
/// as a ``quantum_kernel`` function.
///
/// The data for 2 different molecules are included, and the desired molecule
/// can be selected by changing out the commented variables ``ham_filename``,
/// ``target_energy``, and ``output_fname``.
///
/// Compiles with Intel Quantum SDK version 1.1.3 and later:
/// /<your>/<path>/<to>/intel/quantum-sdk/intel-quantum-compiler -I /<your>/<path>/<to>/ensmallen/include -I /<your>/<path>/<to>/quantum-application-building-blocks/include /<your>/<path>/<to>/quantum-application-building-blocks/examples/vqe_qk.cpp
///
//===----------------------------------------------------------------------===//

// C++ standard library
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Intel(R) Quantum SDK header files
#include <clang/Quantum/quintrinsics.h>
#include <quantum_full_state_simulator_backend.h>

// qabbl/include/
#include <hamiltonian_parsing.h>

// Armadillo https://arma.sourceforge.net/docs.html
//    Conrad Sanderson and Ryan Curtin.
//    Armadillo: An Efficient Framework for Numerical Linear Algebra.
//    International Conference on Computer and Automation Engineering, 2025.
//    Conrad Sanderson and Ryan Curtin.
//    Practical Sparse Matrices in C++ with Hybrid Storage and Template-Based Expression Optimisation.
//    Mathematical and Computational Applications, Vol. 24, No. 3, 2019.
#include <armadillo>

// Ensmallen numerical library https://ensmallen.org/
//   Ryan R. Curtin, Marcus Edel, Rahul Ganesh Prabhu, Suryoday Basak, Zhihao Lou, Conrad Sanderson.
//   The ensmallen library for flexible numerical optimization.
//   Journal of Machine Learning Research, Vol. 22, No. 166, 2021.
#include <ensmallen.hpp>



/////////////////////////////////
// User Inputted constants
//   qubit and circuit parameters for ansatz_hardware_efficient()
const unsigned n_qubits = 4;
const unsigned n_layers = 6;
const unsigned n_params = 3 * n_qubits * n_layers;
//
/*   file where the Hamiltonian for evaluation is saved: */
std::string ham_filename = "hamiltonians/H2/jw_4.ham"; // H2
//std::string ham_filename = "hamiltonians/LiH/jw_4.ham"; // LiH
//
/*   energy of the molecule */
double target_energy = -1.1314597428; // H2
//double target_energy = -7.9836996944; // LiH
//
/*   filename for saved output */
std::string output_fname = "H2_energies.txt"; // H2
//std::string output_fname = "LiH_energies.txt"; // LiH

/////////////////////////////////
// Global Variables
qbit qreg[n_qubits];
double param_reg[n_params];

/////////////////////////////////
// Quantum kernels

quantum_kernel void prep_z_all(){
    for (unsigned i = 0; i < n_qubits; i++) {
        PrepZ(qreg[i]);
    }
}

//  Ansatz appears in
//    Abhinav Kandala, Antonio Mezzacapo, Kristan Temme, Maika Takita, Markus Brink, Jerry M. Chow, and Jay M. Gambetta, 
//    Hardware-efficient variational quantum eigensolver for small molecules and quantum magnets. Nature 549, 242 (2017).
quantum_kernel void ansatz_hardware_efficient(){
    for(unsigned l = 0; l < n_layers; l++){
        // single qubit rotations
        for(unsigned i = 0; i < n_qubits; i++) {
            unsigned param_idx = 3*i + 3*l*n_qubits;
            RY(qreg[i], param_reg[param_idx]);
            RZ(qreg[i], param_reg[param_idx + 1]);
            RY(qreg[i], param_reg[param_idx + 2]);
        }
        
        // entanglement
        for (unsigned i = 0; i < n_qubits - 1; i++)
            CNOT(qreg[i], qreg[i+1]);
    }
}

////////////////////////////////////
// Energy evaluator
// ensmallen's optimizations expect to work on a class method named Evaluate
// that accepts arma::mat as a parameter.
class EnergyEvaluator {
public:
    EnergyEvaluator(const std::vector<HamiltonianTerm>& ham, 
                    iqsdk::FullStateSimulator& device)
        : hamiltonian(ham), iqs_device(device) {}

    double Evaluate(const arma::mat& theta) {
        double energy = run_vqe(iqs_device, theta, hamiltonian);
        
        // Log energy to file
        std::ofstream file(output_fname, std::ios_base::app);
        file << energy << "\n";
        file.close();
        
        return energy;
    }

private:
    const std::vector<HamiltonianTerm>& hamiltonian;
    iqsdk::FullStateSimulator& iqs_device;
    double run_vqe(iqsdk::FullStateSimulator &iqs_device, 
                   const arma::mat &params,
                   const std::vector<HamiltonianTerm>& hamiltonian);
};

// VQE energy calculation
double EnergyEvaluator::run_vqe(
        iqsdk::FullStateSimulator &iqs_device, 
        const arma::mat &params,
        const std::vector<HamiltonianTerm>& hamiltonian) {

    double total_energy = 0.0;

    // set ansatz parameters
    for (unsigned i = 0; i < n_params; i++){
        param_reg[i] = params[i];
    }

    if (iqsdk::QRT_ERROR_SUCCESS != iqs_device.ready()) {
        std::cerr << "IQS initialization failed!" << std::endl;
        return -1;
    }
    // prepare quantum state
    prep_z_all();
    ansatz_hardware_efficient();

    // setup qrefs
    std::vector<std::reference_wrapper<qbit>> qids;
    for (unsigned q = 0; q < n_qubits; q++) {
        qids.push_back(std::ref(qreg[q]));
    }

    // calculate expectation value
    for (const auto& term : hamiltonian) {
        double expectation = iqs_device.getExpectationValue(qids, term.pauli);
        total_energy += term.coeff * expectation;
    }

    return total_energy;
}



/// @cond
int main() {
    // Load Hamiltonian from file
    std::cout << "Loading Hamiltonian from: " << ham_filename << std::endl;
    std::vector<HamiltonianTerm> hamiltonian = load_hamiltonian(ham_filename);
    std::cout << "Successfully loaded Hamiltonian with " << hamiltonian.size() << " terms:\n";
    
    // validate Hamiltonian with n_qubits
    for (const auto& term : hamiltonian) {
        if (term.pauli.length() != n_qubits) {
            std::cerr << "Error: Pauli string length (" << term.pauli.length() 
                      << ") doesn't match n_qubits (" << n_qubits << ")" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    
    // // Print loaded terms
    // std::cout << std::fixed << std::setprecision(16);
    // for (const auto& term : hamiltonian) {
    //     std::cout << term.coeff << " * " << term.pauli << std::endl;
    // }
    // std::cout << std::defaultfloat << std::setprecision(6);
    // std::cout << std::endl;

    // Clear output file
    std::ofstream file(output_fname, std::ios_base::out);
    file.close();

    // Initialize parameters
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.2, 0.2);

    // Setup IQS
    iqsdk::IqsConfig iqs_config(n_qubits, "noiseless");
    iqsdk::FullStateSimulator iqs_device(iqs_config);

    // Create energy evaluator
    EnergyEvaluator evaluator(hamiltonian, iqs_device);

    ////////////////////////////////////////
    // Optimizer Selection

    // SPSA
    ens::SPSA optimizer(0.05, 0.602, 0.1, 0.101, 1000, 1e-8);

    // Genetic Algorithm  
    // ens::CNE optimizer(200, 1000, 1e-12, 50, 0.1);
    
    ////////////////////////////////////////
    // Optimize circuit angles for lowest total energy
    std::cout << "Starting VQE optimization with " << n_params << " parameters...\n";
    std::cout << "Circuit structure: " << n_layers << " layers, " << n_qubits << " qubits\n";

    double optimal_energy = 0.0;
    double best_error = std::numeric_limits<double>::max();
    arma::mat best_params;

    for (unsigned trial = 0; trial < 10; trial++) {
        std::vector<double> init_params(n_params);
        for (unsigned i = 0; i < n_params; i++) {
            init_params[i] = dis(gen);
        }
        
        arma::mat params_trial(init_params);
        double energy = optimizer.Optimize(evaluator, params_trial);  // Generic call
        double error = std::abs(energy - target_energy);
        
        std::cout << "Trial " << trial+1 << ": " << energy 
                  << " Hartree (Error: " << error << ")" << std::endl;
        
        if (error < best_error) {
            best_error = error;
            optimal_energy = energy;
            best_params = params_trial;
        }
    }

    std::cout << "Optimal result: " << optimal_energy 
              << " Hartree (Error: " << best_error << ")" << std::endl;

    // Print results
    std::cout << "\n====== VQE Results ======" << std::endl;
    std::cout << "Optimized Parameters:" << std::endl;
    best_params.print();
    std::cout << "Optimal ground state energy: " << optimal_energy << " Hartree" << std::endl;
    std::cout << "Target ground state energy:  " << target_energy  << " Hartree" << std::endl;
    std::cout << "Error: " << best_error << " Hartree" << std::endl;

    return 0;
}
/// endcond
