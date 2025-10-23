//===--- state_prep.h -------------------------------------------*- C++ -*-===//
//
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
/// \file state_prep.h
/// \brief A quantum expression for preparing a specified state.
///
/// The top-level function is:
///        1. \c prepareState - Takes as input a \c QList and a state
///           (complex vector of size \c 2^N)
///           and returns a \c QExpr that prepares that state.
///
/// Implementation is based on the paper:
///
///   Shende, V. V., S. S. Bullock, and I. L. Markov. "Synthesis of
///   quantum-logic circuits." IEEE Transactions on Computer-Aided Design of
///   Integrated Circuits and Systems 25.6 (2006): 1000-1010.
///
//===----------------------------------------------------------------------===//

#ifndef STATE_PREP_H
#define STATE_PREP_H

// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/quintrinsics.h>
#include <qexpr_utils.h>
#include <qrt_indexing.hpp>

// qabbl/include/
#include <pauli_ops.h> // singleQubitRotation
#include <umatrix.h>

// C++ standard library
#include <cassert>
#include <iostream>
#include <vector>
#include <complex>

/// \var C
/// @brief A type definition for \c std::complex<double>.
using C = std::complex<double>;

const double EULER = 2.71828182845904523536;

using QssIndex = iqsdk::QssIndex;

//////////////////
// Declarations //
//////////////////

/** @addtogroup qstates Quantum States
 * Functions representing quantum states with specific identities and
 * extrapolations.
 * @{
 */

/// @brief Prepare the qubits \c qs in the n-qubit state \c |psi>
///
/// @param qs  A \c QList of size \c N .
/// @param psi A \c 2^N size vector of complex numbers representing a normalized
///            multi-qubit state.
/// @param debug perform consistency checks.
template <int numQubits>
QExpr prepareState(qlist::QList qs, const StateVector<numQubits> &psi,
                   bool debug = false);

/** @} */ // end of qstates


/// @brief Prepares the qubit \c q in the state \f$alpha |0> + beta |1>\f$
///
/// @param q     The \c qbit variable the rotation operates on.
/// @param alpha The amplitude of \f$|0> \f$ in the resulting state.
/// @param beta  The amplitude of \f$|1> \f$ in the resulting state.
/// @param debug perform consistency checks.
QExpr prepareState(qbit &q, C alpha, C beta, bool debug = false);


/// @brief Applies a QMUX of two rotations to qubits \c q0 and \c q1
///        implementing
///        \f$ Rot(P,theta0) \oplus Rot(P,theta1) \f$.
///        Equivalent to
///        `qIf(q0, singleQubitRotation(p,q1,theta0),
///        singleQubitRotation(p,q1,theta1))`.
///
/// @param p             A \c DataList of length 1 equal to "X", "Y", or "Z" to
///                      specify the rotation axis.
/// @param q0            The \c qbit variable to act as the control.
/// @param q1            The \c qbit variable the rotation operates on.
/// @param theta0,theta1 The possible rotation angles
QExpr rotationQMUX(datalist::DataList p, qbit &q0, qbit &q1, double theta0,
                   double theta1);


/// @brief Implements a controlled p-rotation
///        controlled by \c qs
///        where \c qs in state \c |bitstring> yields a rotation by
///        \c theta.
///
/// @tparam numQubits  The number of qubits in the register that the quantum
///                    expression acts on.
/// @param p      A \c DataList of length 1 equal to "X", "Y", or "Z" to specify
///               the rotation axis.
/// @param qs      A \c QList of size \c N of \c qbit variables to act as
///                control qubits.
/// @param q       The \c qbit variable the rotations operate on.
/// @param thetas  A vector of rotation angles of size \c 2^N
/// @param debug   Print additional output for debugging purposes.
///
// See Theorem 8 [ref 1].
template <int numQubits>
PROTECT QExpr rotationQMUX(datalist::DataList p, qlist::QList qs, qbit &q,
                           const std::vector<double> &thetas,
                           bool debug = false);


/// @brief Convert bra-ket notation \f$ alpha|0> + beta|1> \f$
///    to spherical coordinates.
///
/// @returns (r, t, theta, phi) such that
///      \f$ alpha = r * e^{-i t/2} e^{ i phi/2} cos(theta/2) \f$
///      \f$ beta  = r * e^{-i t/2} e^{-i phi/2} sin(theta/2) \f$
///
/// @param alpha The amplitude of \f$|0> \f$ in the input state.
/// @param beta  The amplitude of \f$|1> \f$ in the input state.
/// @param debug if \c true, uses asserts to verify the following holds on the
///              results:
///     \f$ RY(theta) RZ(phi) (alpha |0> + beta |1>) = r e^{-i t / 2} |0> \f$
std::tuple<double, double, double, double>
braketToSpherical(C alpha, C beta, bool debug = false);


////////////////////
// Implementation //
////////////////////

// =======================================================
// Converting between bra-ket and spherical coordinates //
// =======================================================
/// reminder of convenience variable for sqrt(-1) from umatrix.h
//  const C i(0.0, 1.0);
// Want to check that:
//      RY(theta) RZ(phi) (alpha |0> + beta |1>) = r e^{-i t / 2} |0>
void checkBraketToSphericalResults(C alpha, C beta, double r, double t,
                                   double theta, double phi) {
  StateVector<1> state{{alpha, beta}};
  StateVector<1> LHS = RYMatrix(theta) * RZMatrix(phi) * state;
  StateVector<1> RHS{{r * std::exp(-i * t / 2.0), 0.0}};

  if (!isEqual<1>(LHS, RHS)) {
    std::cout << "==ERROR in braketToSpherical==\n\n";
    std::cout << "RY(theta) * RZ(theta) * (alpha,beta) = ( " << LHS[0] << " , "
              << LHS[1] << " )\n";
    std::cout << "r * e^{-i t/2} |0> = ( " << RHS[0] << " , " << RHS[1]
              << " )\n";
    std::cout << "\n";
    std::cout << "In other words, "
              << "(alpha, beta) should be equal to r e^{-i t/2} RZ(-phi) "
                 "RY(-theta) |0>: "
              << r * std::pow(EULER, -i * t / 2.0) * RZMatrix(-phi) *
                     RYMatrix(-theta) * indexToEigen<1>(iqsdk::QssIndex("|0>"))
              << "\n";
    assert(false);
  }
}


std::tuple<double, double, double, double> braketToSpherical(C alpha, C beta,
                                                             bool debug) {
  const double zero_tol = 1e-10; // nearly == 0 special cases

  if (debug) {
    std::cout << "Calling braketToSpherical on alpha = " << alpha
              << ", beta = " << beta << "\n";
  }

  // Case: alpha=beta=0
  if (std::norm(alpha) <= zero_tol && std::norm(beta) <= zero_tol) {
    // r=0; selection of theta, phi, and t should not matter
    double theta = 0.0;
    double phi = 0.0;
    double t = 0.0;
    double r = 0.0;

    if (debug) {
      std::cout << "(case 0) Got r = " << r << ", t = " << t
                << ", theta = " << theta << ", phi = " << phi << "\n";
      checkBraketToSphericalResults(alpha, beta, r, t, theta, phi);
    }
    return std::make_tuple(r, t, theta, phi);
  }

  // Case:  general solution
  // rewrite beta to facilitate
  //   alpha = r * e^{-i t/2} e^{ i phi/2} cos(theta/2)
  //
  //   beta  = r * e^{-i t/2} e^{-i phi/2} (-sin(theta/2))
  //         = r * e^{-i t/2} e^{-i phi/2} e^{i PI} sin(theta/2)
  double r = std::sqrt(std::norm(alpha) + std::norm(beta));
  double theta =
    2.0 * std::atan2(std::sqrt(std::norm(beta)), std::sqrt(std::norm(alpha)));

  double alpha_ang = std::arg(alpha); // == -t/2 + phi/2
  double beta_ang = std::arg(beta);   // == -t/2 - phi/2 + PI
  double phi = M_PI + alpha_ang - beta_ang;
  double t = M_PI - alpha_ang - beta_ang;

  if (debug) {
    std::cout << "(general case) Got r = " << r << ", t = " << t
              << ", theta = " << theta << ", phi = " << phi << "\n";
    checkBraketToSphericalResults(alpha, beta, r, t, theta, phi);
  }
  return std::make_tuple(r, t, theta, phi);
}


// =================================
// Single-qubit state preparation //
// =================================

// Theorem 2 [ref 1]: Given a normalized 1-qubit state |phi>, produce a unitary
//                    circuit U such that
//                    U|phi> = |0>
QExpr disentangle1(qbit &q, C alpha, C beta, bool debug) {
  auto [r, t, theta, phi] = braketToSpherical(alpha, beta, debug);
  return qexpr::_RY(q, theta) * qexpr::_RZ(q, phi);
}


QExpr prepareState(qbit &q, C alpha, C beta, bool debug) {
  return qexpr::_PrepZ(q) + qexpr::invert(disentangle1(q, alpha, beta, debug));
}


// ==============================
// QMUX -- Quantum multiplexer //
// ==============================

template <typename T> void prettyVector(const std::vector<T> &v) {
  std::cout << "\t(\t";
  for (int i = 0; i < v.size(); i++) {
    std::cout << v.at(i) << "\t";
  }
  std::cout << ")\n";
}


PROTECT QExpr rotationQMUX(datalist::DataList p, qbit &q0, qbit &q1,
                           double theta0, double theta1) {
  return singleQubitRotation(p, q1, (theta0 + theta1) / 2.0) +
         qexpr::_CNOT(q0, q1) +
         singleQubitRotation(p, q1, (theta0 - theta1) / 2.0) +
         qexpr::_CNOT(q0, q1);
}


template <>
PROTECT QExpr rotationQMUX<0>(datalist::DataList p, qlist::QList qs, qbit &q,
                              const std::vector<double> &thetas, bool debug) {
  if (debug) {
    std::cout << "Base case: rotationQMUX called with angles ";
    prettyVector(thetas);
    std::cout << "\n";
  }

  assert(qs.size() == 0);
  assert(thetas.size() == 1);
  return singleQubitRotation(p, q, thetas[0]);
}


template <int numQubits>
PROTECT QExpr rotationQMUX(datalist::DataList p, qlist::QList qs, qbit &q,
                           const std::vector<double> &thetas, bool debug) {
  const int newAnglesSize = powerOfTwo(qs.size() - 1);

  if (debug) {
    std::cout << "Recursive case: rotationQMUX called with angles ";
    prettyVector(thetas);
    std::cout << "\n";
  }

  if (numQubits != qs.size() || thetas.size() != powerOfTwo(qs.size())) {
    std::cout << "ERROR: rotationQMUX called with N=" << numQubits << " and "
              << qs.size() << " qubits but only " << thetas.size()
              << " angles, but expects 2^numQubits = " << powerOfTwo(qs.size())
              << " angles.\n";
    assert(0);
  }

  std::vector<double> thetas0; //(newAnglesSize);
  std::vector<double> thetas1; //(newAnglesSize);

  // For the recursive case, the vector thetas (of size 2^numQubits)
  // should be broken up into newthetas0 and newthetas1,
  // both of size 2^(numQubits-1).

  // Use QssIndex patterns to get all indices associated with
  // |0>|phi> and |1>|phi>
  std::vector<int> pattern0(numQubits, QssIndex::Wildcard);
  std::vector<int> pattern1(numQubits, QssIndex::Wildcard);
  assert(pattern0.size() == numQubits);
  assert(pattern1.size() == numQubits);
  pattern0[0] = 0;
  pattern1[0] = 1;
  for (QssIndex idx : QssIndex::patternToIndices(pattern0)) {
    thetas0.push_back(thetas[idx.getIndex()]);
  }
  for (QssIndex idx : QssIndex::patternToIndices(pattern1)) {
    thetas1.push_back(thetas[idx.getIndex()]);
  }

  // Ensure that the pattern-generated vectors have the expected size
  assert(thetas0.size() == newAnglesSize);
  assert(thetas1.size() == newAnglesSize);
  assert(thetas0.size() == thetas1.size());

  std::vector<double> newthetas0(newAnglesSize); // = 0.5 thetas0 + 0.5 thetas1
  std::vector<double> newthetas1(newAnglesSize); // = 0.5 thetas0 - 0.5 thetas1
  for (int i = 0; i < thetas0.size(); i++) {
    assert(i < newthetas0.size() && i < newthetas1.size()); // Bounds check
    newthetas0[i] = 0.5 * (thetas0[i] + thetas1[i]);
    newthetas1[i] = 0.5 * (thetas0[i] - thetas1[i]);
  }

  return
    // recursive case: qs = (qs<<1) + qs[numQubits-1]
    // equivalent to:
    //    qexpr::qIf(qs[numQubits-1],
    //        rotationQMUX<numQubits-1>(p, qs<<1, q, thetas1, debug),
    //        rotationQMUX<numQubits-1>(p, qs<<1, q, thetas0, debug)
    //    );

    rotationQMUX<numQubits - 1>(p, qs << 1, q, newthetas0, debug) +
    qexpr::_CNOT(qs[numQubits - 1], q) +
    rotationQMUX<numQubits - 1>(p, qs << 1, q, newthetas1, debug) +
    qexpr::_CNOT(qs[numQubits - 1], q);
}


/** @addtogroup qstates Quantum States
 * Functions representing quantum states with specific identities and
 * extrapolations.
 * @{
 */

/// @brief     Return a quantum expression that prepares the specified QssIndex.
/// @tparam numQubits  The number of qubits in the register.
/// @param qs  The qubits to be prepared in the specified state.
/// @param idx The state for preparation.
template <int numQubits> PROTECT QExpr prepareQssIndex(qlist::QList qs, QssIndex idx) {
  assert(qs.size() == numQubits);
  assert(idx.getNumQubits() == numQubits);
  double thetas[numQubits];
  for (size_t i = 0; i < numQubits; i++) {
    thetas[i] = idx.bitAt(i) * M_PI;
  }

  return qexpr::map(qexpr::_PrepZ, qs) + qexpr::map2(qexpr::_RX, qs, thetas);
}

/** @} */ // end of qstates


// ============================
// N-qubit state preparation //
// ============================

template <int numQubits>
std::pair<std::vector<double>, std::vector<double>>
partitionAngles(std::vector<double> thetas) {
  std::vector<double> thetas0;
  std::vector<double> thetas1;

  // Use QssIndex patterns to get all indices associated with
  // |0>|phi> and |1>|phi>
  std::vector<int> pattern0(numQubits, QssIndex::Wildcard);
  std::vector<int> pattern1(numQubits, QssIndex::Wildcard);
  pattern0[0] = 0;
  pattern1[0] = 1;
  for (QssIndex idx : QssIndex::patternToIndices(pattern0)) {
    thetas0.push_back(thetas[idx.getIndex()]);
  }
  for (QssIndex idx : QssIndex::patternToIndices(pattern1)) {
    thetas1.push_back(thetas[idx.getIndex()]);
  }
  return std::pair(thetas0, thetas1);
}


// Returns a QExpr U satisfying U |psi> = |0..0>
// Input:
//    numQubits number of qubits satisfying N > 0
//    qs QList of qubits of length numQubits
//    psi state vector of size 2^numQubits
//    debug optional debugging flag
template <int numQubits>
// template recursive case
QExpr disentangleState(qlist::QList qs, const StateVector<numQubits> &psi,
                       bool debug) {
  if (debug) {
    std::cout << "Calling disentangleState with psi =\n" << psi << "\n";
  }
  assert(qs.size() == numQubits);
  assert(psi.rows() == powerOfTwo(numQubits)); // Ensure psi has correct size

  // recursive case
  StateVector<numQubits - 1> next_psi;
  std::vector<double> thetas(next_psi.rows());
  std::vector<double> phis(next_psi.rows());

  // Ensure the vectors have the expected size for rotationQMUX
  assert(next_psi.rows() == powerOfTwo(numQubits - 1));
  assert(thetas.size() == powerOfTwo(numQubits - 1));
  assert(phis.size() == powerOfTwo(numQubits - 1));

  for (unsigned int j = 0; j < next_psi.rows(); j++) {
    // Ensure we don't access out of bounds in the psi vector
    assert(2 * j < psi.rows());
    assert((2 * j) + 1 < psi.rows());

    auto [r_c, t_c, theta_c, phi_c] =
        braketToSpherical(psi(2 * j), psi((2 * j) + 1), debug);

    C psi_entry = r_c * std::exp(-i * t_c / 2.0); // i is complex
    next_psi(j) = psi_entry;
    thetas[j] = theta_c;
    phis[j] = phi_c;
  }

  // disentangleLastQubit<numQubits-1>(qs<<1, qs[numQubits-1], thetas, phis)
  return rotationQMUX<numQubits - 1>("Z", qs >> 1, qs[0], phis, debug) +
         rotationQMUX<numQubits - 1>("Y", qs >> 1, qs[0], thetas, debug) +
         disentangleState<numQubits - 1>(qs >> 1, next_psi, debug);
}


template <>
// template base case
QExpr disentangleState<1>(qlist::QList qs, const StateVector<1> &psi,
                          bool debug) {
  assert(qs.size() == 1);
  return disentangle1(qs[0], psi[0], psi[1], debug);
}


template <int numQubits>
QExpr prepareState(qlist::QList qs, const StateVector<numQubits> &psi,
                   bool debug) {
  assert(qs.size() == numQubits);
  return qexpr::map(qexpr::_PrepZ, qs) +
         qexpr::invert(disentangleState<numQubits>(qs, psi, debug));
}

#endif
