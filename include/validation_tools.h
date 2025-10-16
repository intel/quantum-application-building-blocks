//===--- validation_tools.h -------------------------------------*- C++ -*-===//
//
//===----------------------------------------------------------------------===//
//
// Copyright 2025 Intel Corporation.
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
/// \file validation_tools.h
/// \brief A library of tools to validate the behavior of quantum expressions.
///
/// A set of functions to compare the output from the Intel Quantum Simulator
/// to equivalent analytic expressions to validate the behavior of quantum
/// algorithms written as QExpr type.
///
//===----------------------------------------------------------------------===//

#ifndef IQSDK_VALID_H
#define IQSDK_VALID_H

using namespace qlist;
using namespace iqsdk;
using namespace qexpr;


//////////////////
// Declarations //
//////////////////

/**
 * @defgroup validation Validation tools
 * A library of functions for comparing simulated qubit results to analytic
 * expressions.
 * @{
 */

/// @brief             Returns a bool that signals if the state produced by
///                    executing e is equivalent to psi.
/// @warning           This function does not perform initialization unless the
///                    steps of e explicitly performs initialization.
/// @tparam numQubits  The number of qubits in the expected StateVector and
///                    operated on by the QExpr to be checked.
/// @param iqs_device  A ready FullStateSimulator class to simulate the quantum
///                    expressions.
/// @param qs          The QList of qubits to run the quantum expression on.
/// @param e           The QExpr to be checked.
/// @param psi         The expected state of the qubit register after applying
///                    the QExpr e.
/// @param debug       Perform internal checks and print amplitudes and
///                    probabilities for debugging purposes
template <int numQubits>
__attribute__((always_inline)) bool
implementsOutputState(FullStateSimulator &iqs_device, QList qs, QExpr e,
                      const StateVector<numQubits> &psi, bool debug = false) {
  eval_hold(e);
  // eval_hold( qexpr::printQuantumLogic(e));

  // get the output state and compare to psi
  auto qsRefs = to_ref_wrappers(qs);
  auto amps = iqs_device.getAmplitudes(qsRefs, {});
  StateVector<numQubits> output1 = qssMapToEigen<numQubits>(amps);
  bool result = isEqualUpToPhase<numQubits>(output1, psi);

  if (debug && !result) {
    std::cout << "Failed implementsOutputState with \n";

    std::cout << "Target state:\n" << psi << "\n";
    auto expected_amps = eigenToQssMap<numQubits>(psi);
    FullStateSimulator::displayAmplitudes(expected_amps);
    std::cout << "Simulator amplitudes:\n";
    FullStateSimulator::displayAmplitudes(amps);

    std::cout << "Target probabilities:\n";
    auto expected_probs =
        FullStateSimulator::amplitudesToProbabilities(expected_amps);
    FullStateSimulator::displayProbabilities(expected_probs);
    std::cout << "Simulator probabilities:\n";
    auto actual_probs = FullStateSimulator::amplitudesToProbabilities(amps);
    FullStateSimulator::displayProbabilities(actual_probs);
  }

  return result;
}


/// @brief             Returns a bool that signals if the state produced by
///                    preparing psi followed by executing e is equivalent
///                    applying the matrix U to psi.
/// @tparam numQubits  The number of qubits in the expected StateVector and
///                    operated on by the QExpr to be checked.
/// @param iqs_device  A ready FullStateSimulator class to simulate the quantum
///                    expressions.
/// @param qs          The QList of qubits to run the quantum expression on.
/// @param psi         The initial state of the qubit register to prepare before
///                    applying the QExpr e.
/// @param e           The QExpr to be checked.
/// @param U           Matrix representation of the QExpr e in UMatrix form.
/// @param debug       Perform internal checks and print amplitudes and
///                    probabilities for debugging purposes
template <int numQubits>
__attribute__((always_inline)) bool
implementsUMatrixOn(FullStateSimulator &iqs_device, QList qs,
                    const StateVector<numQubits> &psi, QExpr e,
                    const UMatrix<numQubits> &U, bool debug = false) {
  eval_hold(prepareState<numQubits>(qs, psi, false));

  auto qsRefs = to_ref_wrappers(qs);
  auto probs = iqs_device.getProbabilities(qsRefs, {});
  if (debug) {
    std::cout << "Prepared state: ";
    FullStateSimulator::displayProbabilities(probs);
  }

  StateVector<numQubits> output_psi = U * psi;

  bool result = implementsOutputState<numQubits>(iqs_device, qs, e, output_psi,
                                                 debug);
  if (debug && !result) {
    std::cout << "Failed implementsUMatrixOn with \n";
    std::cout << "Input state:\n" << psi << "\n";
    std::cout << "U:\n" << U << "\n";
  }
  return result;
}


/// @brief             Returns a bool that signals if the state produced by
///                    preparing basis state idx followed by executing e is
///                    equivalent to applying the matrix U to idx.
/// @tparam numQubits  The number of qubits in the expected StateVector and
///                    operated on by the QExpr to be checked.
/// @param iqs_device  A ready FullStateSimulator class to simulate the quantum
///                    expressions.
/// @param qs          The QList of qubits to run the quantum expression on.
/// @param idx         The basis state of the qubit register to prepare before
///                    applying the QExpr e.
/// @param e           The QExpr to be checked.
/// @param U           Matrix representation of the QExpr e in UMatrix form.
/// @param debug       Perform internal checks and print amplitudes and
///                    probabilities for debugging purposes
template <int numQubits>
__attribute__((always_inline)) bool
implementsUMatrixOn(FullStateSimulator &iqs_device, QList qs,
                    QssIndex idx, QExpr e, const UMatrix<numQubits> &U,
                    bool debug = false) {
  eval_hold(prepareQssIndex<numQubits>(qs, idx));
  StateVector<numQubits> psi = indexToEigen<numQubits>(idx);

  auto qsRefs = to_ref_wrappers(qs);
  auto probs = iqs_device.getProbabilities(qsRefs, {});
  if (debug) {
    std::cout << "Prepared state: ";
    FullStateSimulator::displayProbabilities(probs);
  }

  StateVector<numQubits> output_psi = U * psi;

  bool result = implementsOutputState<numQubits>(iqs_device, qs, e, output_psi,
                                                 debug);
  if (debug && !result) {
    std::cout << "Failed implementsUMatrixOn with \n";
    std::cout << "Input state: " << idx << "\n";
    std::cout << "U:\n" << U << "\n";
  }
  return result;
}


/// @brief             Returns a bool that signals if the state produced by
///                    preparing a pseudo-random state |psi> and executing
///                    QExpr e is equivalent to the state produced by applying
///                    U to |psi>
/// @tparam numQubits  The number of qubits in the expected StateVector and
///                    operated on by the QExpr to be checked.
/// @param iqs_device  A ready FullStateSimulator class to simulate the quantum
///                    expressions.
/// @param generator   A \c std::mt19937 random number generator.
/// @param num_tries   The number of pseudo-random states to test on.
/// @param qs          The QList of qubits to run the quantum expression on.
/// @param e           The QExpr to be checked.
/// @param U           Expected matrix representation of the \c QExpr e.
/// @param debug       Perform internal checks and print amplitudes and
///                    probabilities for debugging purposes
template <int numQubits>
__attribute__((always_inline)) bool
implementsUMatrixHaar(FullStateSimulator &iqs_device,
                      std::mt19937 &generator, int num_tries, QList qs, QExpr e,
                      const UMatrix<numQubits> &U, bool debug = false) {
  for (int i = 0; i < num_tries; i++) {
    StateVector<numQubits> psi = haar_state_distribution<numQubits>(generator);
    if (!implementsUMatrixOn<numQubits>(iqs_device, qs, psi, e, U, debug)) {
      std::cout << "FAILED on iteration " << i << "\n";
      return false;
    }
  }
  if (debug) {
    std::cout << "Successfully tested implementsUMatrix on " << num_tries
              << " random state(s)\n";
  }
  return true;
}


/// @brief             Returns a bool that signals if the state produced by
///                    executing QExpr e is equivalent to applying U for all
///                    basis states in the qubit register.
/// @tparam numQubits  The number of qubits in the expected StateVector and
///                    operated on by the QExpr to be checked.
/// @param iqs_device  A ready FullStateSimulator class to simulate the quantum
///                    expressions.
/// @param qs          The QList of qubits to run the quantum expression on.
/// @param e           The QExpr to be checked.
/// @param U           Expected matrix representation of the \c QExpr e.
/// @param debug       Perform internal checks and print amplitudes and
///                    probabilities for debugging purposes
template <int numQubits>
__attribute__((always_inline)) bool
implementsUMatrix(FullStateSimulator &iqs_device, QList qs, QExpr e,
                  const UMatrix<numQubits> &U, bool debug = false) {
  for (unsigned int i = 0; i < powerOfTwo(numQubits); i++) {
    QssIndex idx(numQubits, i);
    StateVector<numQubits> phi = indexToEigen<numQubits>(idx);
    bool b = implementsUMatrixOn<numQubits>(iqs_device, qs, idx, e, U, debug);

    if (!b) {
      std::cout << "Failed implementsUMatrix\n";
      return false;
    }
  }
  return true;
}

/** @} */ // end of implementsNNN


/// @brief             Prints the probabilities and amplitudes of the qubit
///                    register.
/// @param iqs_device  A ready FullStateSimulator object to draw results from.
/// @param qsRefs      A set of references to the qubits that form the qubit
///                    register or state-space of interest.
void displayCurrentProbs(FullStateSimulator &iqs_device,
                         std::vector<std::reference_wrapper<qbit>> &qsRefs) {

  auto probs = iqs_device.getProbabilities(qsRefs, {});
  FullStateSimulator::displayProbabilities(probs);

  auto amps = iqs_device.getAmplitudes(qsRefs, {});
  FullStateSimulator::displayAmplitudes(amps);
}


/// @brief             Throw an error if the current state of the simulator is
///                    not equal to phi (up to a phase)
/// @tparam numQubits  The number of qubits in the expected StateVector and
///                    operated on by the QExpr to be checked.
/// @param iqs_device  A ready FullStateSimulator object with a state to check.
/// @param qsRefs      A set of references to the qubits that form the qubit
///                    register or state-space of interest.
/// @param phi         A StateVector with the expected state of the iqs_device.
/// @param display     Select to display the amplitudes and probabilities.
template <int numQubits>
void assertState(FullStateSimulator &iqs_device,
                 std::vector<std::reference_wrapper<qbit>> &qsRefs,
                 const StateVector<numQubits> &phi, bool display = false) {
  auto amps = iqs_device.getAmplitudes(qsRefs, {});
  auto phi_actual = qssMapToEigen<numQubits>(amps);

  if (!isEqualUpToPhase<numQubits>(phi, phi_actual)) {
    auto probs_actual = FullStateSimulator::amplitudesToProbabilities(amps);

    auto amps_target = eigenToQssMap<numQubits>(phi);
    auto probs_target =
        FullStateSimulator::amplitudesToProbabilities(amps_target);

    std::cout << "ERROR\n"
              << "Expected probabilities:\n";
    FullStateSimulator::displayProbabilities(probs_target);
    std::cout << "Simulator probabilities:\n";
    FullStateSimulator::displayProbabilities(probs_actual);
    std::cout << "\nExpected amplitudes:\n";
    FullStateSimulator::displayAmplitudes(amps_target);
    std::cout << "Simulator amplitudes:";
    FullStateSimulator::displayAmplitudes(amps);
    assert(false);
  }
  if (display) {
    auto amps_target = eigenToQssMap<numQubits>(phi);

    std::cout << "==results of assertState==\n";
    std::cout << "Expected state:\n";
    FullStateSimulator::displayAmplitudes(amps_target);

    std::cout << "Simulator state:\n";
    FullStateSimulator::displayAmplitudes(amps);

    auto probs_target =
        FullStateSimulator::amplitudesToProbabilities(amps_target);
    std::cout << "Expected and target state and amplitudes both produce the "
                 "expected probabilities:\n";
    FullStateSimulator::displayProbabilities(probs_target);
    std::cout << "\n";
  }
}


#endif // IQSDK_VALID_H
