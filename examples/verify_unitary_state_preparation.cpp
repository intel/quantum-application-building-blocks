//===----------------------------------------------------------------------===//
//
// Copyright (C) 2025 Intel Corporation.
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
/// \file verify_unitary_state_preparation.cpp
/// \brief This source code demonstrates preparing qubits into a state.
///
/// Each state prepared on in the qubit simulator is given by the quantum
/// amplitudes. The program checks for the correctness
/// of several hard-coded states for one to three qubits, and then checks the
/// correctness of randomly generated three and four qubit states.
///
/// This example requires features from c++-20 and it needs to be compiled with
/// the appropriate flag given to IQSDK v1.1.1 with the -f prefix, i.e.
/// -f--std=c++20
///
/// e.g.
/// /<path to intel quantum sdk>/intel-quantum-compiler -I<path to repository>/quantum-application-building-blocks/include -I<path to repository>/eigen -f--std=c++20 -o<path to output> <path to repository>/examples/validate_unitary_state_preparation.cpp
/// where the < > represent paths that vary from user to user.
//===----------------------------------------------------------------------===//
// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/quintrinsics.h>
#include <qexpr_utils.h>
#include <qrt_indexing.hpp>
#include <quantum_full_state_simulator_backend.h>

// qabbl/include/
#include <haar_rand.h>
#include <state_prep.h>

// C++ standard library
#include <cassert>
#include <iostream>
#include <iterator>
#include <vector>

/// @brief Test that \c prepareState() has the same probability distribution as
///        \c state produces using an inputted state.
///
/// @tparam qubit_reg_size The number of qubits in the state.
/// @param iqs_device     The simulator object to evaluate the \c Qexpr on.
/// @param state          A vector of \c 2^N complex amplitudes that describe an
///                       \c N qubit state.
/// @param debug          Also produce the consistency check included in
///                       \c prepareState.
///
/// @warning The algorithm in state_prep.h will normalize the input vector,
///          or rescale the elements so that the total length == 1.
///          however, the check implemented in this program does not normalize
///          the input vector, and so expects all input vectors to be normalized
template <int qubit_reg_size>
void testPrepareState(iqsdk::FullStateSimulator &iqs_device,
                      std::vector<C> &state, bool debug = true) {
  assert(state.size() == std::pow(2, qubit_reg_size));

  qbit listable(qs, qubit_reg_size);
  qexpr::eval_hold(prepareState(qs, state, debug));

  std::cout << "== Square of input amplitudes ==" << std::endl;
  auto input_amps = iqsdk::qssVectorToMap(state, qubit_reg_size);
  auto input_probs =
      iqsdk::FullStateSimulator::amplitudesToProbabilities(input_amps);
  iqsdk::FullStateSimulator::displayProbabilities(input_probs);

  std::cout << "== Output probabilities ==" << std::endl;
  auto qubitRefs = to_ref_wrappers(qs);
  auto probs = iqs_device.getProbabilities(qubitRefs, {});
  iqsdk::FullStateSimulator::displayProbabilities(probs);

  // Check equality of vectors
  const double tolerance = 0.1;
  std::vector<double> output_state = qssMapToVector(probs, qubit_reg_size, 0.0);
  assert(state.size() == output_state.size());
  for (size_t idx = 0; idx < state.size(); idx++) {
    assert(std::abs(std::norm(state[idx]) - output_state[idx]) <= tolerance);
  }
}

/// @brief Test that \c prepareState() has the same probability distribution as
///        \c state produces using a randomly generated state.
///
/// @tparam qubit_reg_size The number of qubits in the state.
/// @param iqs_device     The simulator object to evaluate the \c Qexpr on.
///
/// @warning The algorithm in state_prep.h will normalize the input vector,
///          or rescale the elements so that the total length == 1.
///          however, the check implemented in this program does not normalize
///          the input vector, and so expects all input vectors to be normalized
template <int qubit_reg_size>
void testRandomState(iqsdk::FullStateSimulator &iqs_device) {

  // Initialize a vector representing |00...>
  int n_elements = std::pow(2, qubit_reg_size);
  // initialize as all 0's to avoid n memory allocations through push_back
  std::vector<C> computational_zero(n_elements, C(0.0, 0.0));
  computational_zero.at(0) = C(1.0, 0.0);

  // b. Prepare a Haar unitary matrix to apply to |00...>

  // an engine that produces a sequence of pseudo-random values.
  // std::mt19937 engine(1234);  // fixed seed for deterministic sequence
  std::random_device rd;     // random seed
  std::mt19937 engine(rd()); // non-deterministic sequence

  UMatrix matrix_rep = sampleHaarDistributedUnitary(n_elements, engine);

  auto state_to_prep = matrixMultiply(matrix_rep, computational_zero);
  testPrepareState<qubit_reg_size>(iqs_device, state_to_prep, false);
}

/// @cond
int main() {
  const int max_qubits = 4; // problem size:  # of qubits

  iqsdk::IqsConfig iqs_config;
  iqs_config.num_qubits = max_qubits;
  iqsdk::FullStateSimulator iqs_device(iqs_config);
  iqsdk::QRT_ERROR_T status = iqs_device.ready();
  assert(status == iqsdk::QRT_ERROR_SUCCESS);

  // Unit tests
  std::vector<C> state_to_prep;
  double sqrt2 = std::sqrt(2.0);

  // 1-qubit unit tests
  // test 1.1: random state
  testRandomState<1>(iqs_device);
  // test 1.2: |0>
  state_to_prep = {1.0, 0.0};
  testPrepareState<1>(iqs_device, state_to_prep);
  // test 1.3: |0> + |1>
  state_to_prep = {1.0 / sqrt2, 1.0 / sqrt2};
  testPrepareState<1>(iqs_device, state_to_prep);

  // 2-qubit unit tests
  // test 2.1: random state
  testRandomState<2>(iqs_device);
  // test 2.2: |00>
  state_to_prep = {1.0, 0.0, 0.0, 0.0};
  testPrepareState<2>(iqs_device, state_to_prep);
  // test 2.3: |01>
  state_to_prep = {0.0, 0.0, 1.0, 0.0};
  testPrepareState<2>(iqs_device, state_to_prep);
  // test 2.4: |11>
  state_to_prep = {0.0, 0.0, 0.0, 1.0};
  testPrepareState<2>(iqs_device, state_to_prep);
  // test 2.5: |00> + |11>
  state_to_prep = {1.0 / sqrt2, 0.0, 0.0, 1.0 / sqrt2};
  testPrepareState<2>(iqs_device, state_to_prep);
  // test 2.6: |00> + |01> + |10> + |11>
  state_to_prep = {1.0 / (2.0), 1.0 / (2.0), 1.0 / (2.0), 1.0 / (2.0)};
  testPrepareState<2>(iqs_device, state_to_prep);
  // test 2.7: |01> + |10>
  state_to_prep = {0.0, 1.0 / sqrt2, 1.0 / sqrt2, 0.0};
  testPrepareState<2>(iqs_device, state_to_prep);

  // 3-qubit tests
  // test 3.1: random state
  testRandomState<3>(iqs_device);
  // test 3.2: |000> + |111>
  iqsdk::QssMap<C> state_to_prep_qss;
  state_to_prep_qss[iqsdk::QssIndex("|000>")] = 1.0 / sqrt2;
  state_to_prep_qss[iqsdk::QssIndex("|111>")] = 1.0 / sqrt2;
  state_to_prep = qssMapToVector(state_to_prep_qss, 3, {0.0, 0.0});
  testPrepareState<3>(iqs_device, state_to_prep);
  // test 3.3: |000> + |010> + |011> + |111>
  // iqsdk::QssMap<C> state_to_prep_qss;
  state_to_prep_qss.clear();
  state_to_prep_qss[iqsdk::QssIndex("|000>")] = 1.0 / 2.0;
  state_to_prep_qss[iqsdk::QssIndex("|111>")] = 1.0 / 2.0;
  state_to_prep_qss[iqsdk::QssIndex("|010>")] = 1.0 / 2.0;
  state_to_prep_qss[iqsdk::QssIndex("|011>")] = 1.0 / 2.0;
  state_to_prep = qssMapToVector(state_to_prep_qss, 3, {0.0, 0.0});
  testPrepareState<3>(iqs_device, state_to_prep);

  // 4-qubit tests
  // test 4.1:  random state
  testRandomState<4>(iqs_device);

  return 0;
}
/// @endcond
