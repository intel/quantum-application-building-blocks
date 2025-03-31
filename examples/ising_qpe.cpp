//===----------------------------------------------------------------------===//
//
// Copyright (C) 2024 Intel Corporation.
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
/// \file ising_qpe.cpp
/// \brief Demonstrate the quantum phase estimation qpe() .
///
/// Quantum phase estimation
/// (QPE) is a subroutine which takes in any unitary U, and approximates
/// the eigenvalue into a separate register for a prepared eigenstate.
/// In this file, QPE is demonstrated by computing the energy of a quantum system
/// using both a classical and a quantum algorithm and comparing the two results.
///
/// The Quantum Fourier Transform (QFT) is a major subroutine used by QPE.
/// We use QFT without explanation. See "qexpr_qft.cpp" for a demonstration
/// focusing on QFT.
///
/// See "qalgorithm.h" for detailed descriptions.
///
//===----------------------------------------------------------------------===//
// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/quintrinsics.h>
#include <qexpr_utils.h>
#include <quantum_full_state_simulator_backend.h>

// qabbl/include/
#include <qalgorithm.h>

// C++ standard library
#include <bit>
#include <cassert>
#include <iostream>
#include <vector>

/*
Test QPE with a 1D Random-bond Ising model.

Build the random bond instance into the dynamic angles by passing a specific
random-bond instance as an integer which encodes the sign of the Ising terms.
*/

// This function represents a general ZZ Ising term
QExpr ZZPhase(qbit &q1, qbit &q2, double angle, double sign) {
  return qexpr::_CNOT(q1, q2) + qexpr::_RZ(q2, sign * angle) +
         ~qexpr::_CNOT(q1,
                       q2) // NOTE: we need to invert to cancel a residual phase
         + qexpr::global_phase(-angle);
}

// Recursive helper function for generating each Ising term in a 1D chain
QExpr OneDRBIMHelper(unsigned long RBinstance, double angle, qbit &base,
                     qlist::QList chain) {
  // extract a bit of the instance and translate it to a sign
  double sign = (1. - 2. * (double)(RBinstance & 1));
  return qexpr::cIf(
      chain.size() > 0,
      // represents the ZZ term for the begining of the chain and the passed
      // qubit "base".
      ZZPhase(base, chain[0], angle, sign)
          // recurse by passing first qubit as "base" and the slice of "chain".
          + OneDRBIMHelper(RBinstance >> 1, angle, chain[0], chain >> 1),
      qexpr::identity());
}

QExpr OneDRBIM(unsigned long RBInstance, qlist::QList chain, double norm) {
  // need to scale the angle by 1/"norm" and 2 Pi to insure the energy fits
  // inside the phase register.
  double scale = 2 * M_PI / norm;
  // pass the last qubit of the chain as the "base" of the helper to insure
  //  the Ising model forms a cycle.
  return OneDRBIMHelper(RBInstance, scale, chain[chain.size() - 1], chain);
}

// prep classical state at runtime using an integer encoding
QExpr prepTo(qlist::QList qs, unsigned to) {
  unsigned sz = qs.size();
  return qexpr::cIf(sz > 0,
                    prepTo(qs >> 1, to >> 1) *
                        qexpr::_RX(qs[0], M_PI * (to & 1)) *
                        qexpr::_PrepZ(qs[0]),
                    qexpr::identity());
}

// This is a purely classical function which translates a cbit array into a
// "big-endian" integer which is passed by reference.
// we will write the functon with FLEQ so that it can be bound inside the
// "measureTo" QExpr function below.
PROTECT QExpr writeTo(cbit *c, unsigned &out, unsigned shft) {
  out += (*c) << shft;
  return qexpr::cIf(
      shft > 0,
      writeTo(c + 1, out, shft - 1) +
          qexpr::identity(), // NOTE: this is needed to avoid a known bug.
      qexpr::identity());
}

PROTECT QExpr measureTo(qlist::QList reg, unsigned &out) {
  cbit c[reg.size()];
  return qexpr::map(qexpr::_MeasZ, reg, c) << writeTo(c, out, reg.size() - 1);
}

// This will calculate the energy classically from the state, instance and size
// of the chain.
//   To demonstrate, imagine we have the two bit strings, the "state" and the RB
//   "instance":
//
//   state    : 0 1 0 0 1 1 0
//   instance : 1 1 0 1 0 0 1
//
//   the Ising terms are NN parity which we can calculate via a "circular" shift
//   to the right and a bit-wise XOR:
//
//   state    : 0 1 0 0 1 0 0
//   shift    : 0 0 1 0 0 1 0
//            ^ -------------
//   ising    : 0 1 1 0 1 1 0
//
//   Finally we XOR this with the RB instance and the energy is the Ham wieght
//   of the result: ising    : 0 1 1 0 1 1 0 instance : 1 1 0 1 0 0 1
//            ^ -------------
//   energy = | 1 0 1 1 1 1 1 | = 6
unsigned calcEclassical(unsigned long RBInstance, unsigned long state,
                        unsigned sz) {
  unsigned out = 0;
  // right curcular shift
  //  main shift by 1
  unsigned long shift = (state << 1);
  // wrap around
  shift = shift ^ ((state >> (sz - 1)) & 1);
  // XOR for the ising terns
  unsigned long ising = state ^ shift;
  // final energy state
  unsigned long energy = ising ^ RBInstance;
  return std::__popcount(energy & ((1 << sz) - 1));
}

/// @cond
int main() {
  iqsdk::IqsConfig iqs_config;
  iqsdk::FullStateSimulator iqs_device(iqs_config);
  iqsdk::QRT_ERROR_T status = iqs_device.ready();
  assert(status == iqsdk::QRT_ERROR_SUCCESS);

  // Sets the number of try to compare the QPE calculation to the classical
  // calculation.
  unsigned attempts = 20;

  const unsigned chain_size = 15; // requires at least 4 bits to measures
  const unsigned phase_size = 4;

  // state qubits
  qbit listable(chain, chain_size);

  // phase register qubits
  qbit listable(phase, phase_size);

  unsigned success = 0;
  for (int a = 0; a < attempts; a++) {
    // generate random state
    unsigned state = rand() & ((1 << chain_size) - 1);
    // generate random RB instance
    unsigned instance = rand() & ((1 << chain_size) - 1);

    // Use classical method to calculate the energy
    unsigned e_c = calcEclassical(instance, state, chain_size);

    // use QPE to calculate the energy.
    unsigned e_q = 0;
    qexpr::eval_hold(
        prepTo(chain, state) +
        qpe(phase, OneDRBIM(instance, chain, (1 << (phase_size)))) +
        measureTo(phase, e_q));

    // compare classical and quantum energies and report the results.
    if (e_q == e_c) {
      success++;
      std::cout << "QPE SUCCEEDED for instance " << instance << " and state "
                << state << ".\n";
    } else {
      std::cout << "QPE FAILED for instance " << instance << " and state "
                << state << ".\n";
      std::cout << "   "
                << "QPE: " << e_q << " vs classical: " << e_c << "\n";
    }
  }

  std::cout << "QPE was correct " << (double)(success) / (double)(attempts)*100
            << "% of the time.\n";
}
/// @endcond
