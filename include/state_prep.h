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
/// \brief A library for producing quantum kernel expressions that take as input
///        states.
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

// C++ standard library
#include <cassert>
#include <iostream>
#include <vector>

/// \var C
/// @brief A type definition for \c std::complex<double>.
using C = std::complex<double>;

/// \var UMatrix
/// @brief A type definition for building matrices from a \c std::vector of
///        \c std::vector representing the rows. Each elements of the row
///        \c std::vector is a \c std::complex<double>.
using UMatrix = std::vector<std::vector<C>>;

//////////////////
// Declarations //
//////////////////

/// @brief Prepare the qubits \c qs in the n-qubit state \c |psi>
///
/// @param qs  A \c QList of size \c N .
/// @param psi A \c 2^N size vector of complex numbers representing a normalized
///            multi-qubit state.
/// @param debug perform consistency checks.
QExpr prepareState(qlist::QList qs, std::vector<C> psi, bool debug = false);

/// @brief Prepares the qubit \c q in the state \f$alpha |0> + beta |1>\f$
///
/// @param q     The \c qbit variable the rotation operates on.
/// @param alpha The amplitude of \f$|0> \f$ in the resulting state.
/// @param beta  The amplitude of \f$|1> \f$ in the resulting state.
/// @param debug perform consistency checks.
QExpr prepareState(qbit &q, C alpha, C beta, bool debug = false);

/// @brief Applies a single-qubit rotation around the Pauli axis specified
///        by \c p.
///
/// @param p     A \c DataList of length 1 equal to "X", "Y", or "Z" to specify
///              the rotation axis.
/// @param q     The \c qbit variable the rotation operates on.
/// @param theta The rotation angle about the specified axis.
QExpr singleQubitRotation(datalist::DataList p, qbit &q, double theta);

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
/// @param p      A \c DataList of length 1 equal to "X", "Y", or "Z" to specify
///               the rotation axis.
/// @param qs      A \c QList of size \c N of \c qbit variables to act as
///                control qubits.
/// @param q       The \c qbit variable the rotations operate on.
/// @param thetas  A vector of rotation angles of size \c 2^N
///
/// See Theorem 8 [ref 1].
PROTECT QExpr rotationQMUX(datalist::DataList p, qlist::QList qs, qbit &q,
                           std::vector<double> &thetas);

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

/// @brief Return the result of applying the matrix to the input vector
std::vector<C> matrixMultiply(const UMatrix &matrix,
                              const std::vector<C> &vector);

/// @brief Return the product of two matrices.
UMatrix matrixMultiply(const UMatrix &U, const UMatrix &V);

/// @brief Convenience variable for sqrt(-1)
const C i(0.0, 1.0);

////////////////////
// Implementation //
////////////////////

// =========================================
// Complex number and matrix calculations ==
// =========================================

/// Check that two complex numbers are equal
void assert_equal(C alpha, C beta, double tolerance = 0.01) {
  assert(std::abs(std::real(alpha) - std::real(beta)) < tolerance);
  assert(std::abs(std::imag(alpha) - std::imag(beta)) < tolerance);
}

std::vector<C> matrixMultiply(const UMatrix &matrix,
                              const std::vector<C> &vector) {
  std::vector<C> result_vector;
  int dimension = matrix.size();
  assert(matrix.at(0).size() == dimension); // matrix is square
  for (int row = 0; row < dimension; row++) {
    C sum(0.0, 0.0);
    for (int j = 0; j < dimension; j++) {
      sum += matrix.at(row).at(j) * vector.at(j);
    }
    result_vector.push_back(sum);
  }
  return result_vector;
}

UMatrix RYMatrix(double theta) {
  UMatrix RY = {{std::cos(theta / 2.0), std::sin(theta / 2.0)},
                {-std::sin(theta / 2.0), std::cos(theta / 2.0)}};
  return RY;
}
UMatrix RZMatrix(double theta) {
  UMatrix RZ = {{std::exp(-i * theta / 2.0), 0},
                {0, std::exp(i * theta / 2.0)}};
  return RZ;
}

// =======================================================
// Converting between bra-ket and spherical coordinates //
// =======================================================

// Want to check that:
//      RY(theta) RZ(phi) (alpha |0> + beta |1>) = r e^{-i t / 2} |0>
void checkBraketToSphericalResults(C alpha, C beta, double r, double t,
                                   double theta, double phi) {

  std::vector<C> state = {alpha, beta};
  std::vector<C> LHS =
      matrixMultiply(RYMatrix(theta), matrixMultiply(RZMatrix(phi), state));
  std::vector<C> RHS = {r * std::exp(-i * t / 2.0), 0.0};
  std::cout << "LHS: ( " << LHS[0] << " , " << LHS[1] << " )\n";
  std::cout << "RHS: ( " << RHS[0] << " , " << RHS[1] << " )\n";

  assert(LHS.size() == 2);
  assert_equal(LHS[0], RHS[0]);
  assert_equal(LHS[1], RHS[1]);
}

std::tuple<double, double, double, double> braketToSpherical(C alpha, C beta,
                                                             bool debug) {
  const double zero_tol = 1e-30; // nearly == 0 special cases

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

  // Case:  amplitude beta of |1> is zero
  else if (std::norm(beta) <= zero_tol) {

    double theta = 0.0;

    // theta=0 reduces the constraints to the following:
    //      alpha = r * e^{-i t/2} e^{i phi/2}
    //      beta  = 0

    double r = std::sqrt(std::norm(alpha));
    C cphi = -2.0 * i * std::log(alpha / r);
    if (debug) {
      assert(std::abs(std::imag(cphi)) <= zero_tol);
    }
    double phi = std::real(cphi);
    double t = 0.0;

    if (debug) {
      std::cout << "(case 1) Got r = " << r << ", t = " << t
                << ", theta = " << theta << ", phi = " << phi << "\n";
      checkBraketToSphericalResults(alpha, beta, r, t, theta, phi);
    }
    return std::make_tuple(r, t, theta, phi);
  }

  // Case: amplitude alpha of |0> is zero
  else if (std::norm(alpha) <= zero_tol) {

    double theta = M_PI;

    // theta=pi reduces the constraints to the following:
    //    alpha = 0
    //    beta  = r * e^{i t/2} e^{i phi/2}

    double r = std::sqrt(std::norm(beta));
    C cphi = 2.0 * i * std::log(beta / r);
    if (debug) {
      assert(std::abs(std::imag(cphi)) <= zero_tol);
    }
    double phi = std::real(cphi);
    double t = 0.0;

    if (debug) {
      std::cout << "(case 2) Got r = " << r << ", t = " << t
                << ", theta = " << theta << ", phi = " << phi << "\n";
      checkBraketToSphericalResults(alpha, beta, r, t, theta, phi);
    }
    return std::make_tuple(r, t, theta, phi);
  }

  // Case:  general solution
  // all Bloch angles are computed with ratios of alpha and beta so that r is
  // canceled, i.e., the accuracy of the angles doesn't depend on the input
  // satisfying r = 1
  double r = std::sqrt(std::norm(alpha) + std::norm(beta));
  double theta = 2.0 * std::acos(std::abs(alpha) / r);
  double alpha_ang = std::arg(alpha);
  double beta_ang = std::arg(beta);
  double phi = beta_ang - alpha_ang;
  double t = alpha_ang + beta_ang;

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

QExpr rotationQMUX(datalist::DataList p, qbit &q0, qbit &q1, double theta0,
                   double theta1) {
  return singleQubitRotation(p, q1, (theta0 + theta1) / 2.0) +
         qexpr::_CNOT(q0, q1) +
         singleQubitRotation(p, q1, (theta0 - theta1) / 2.0) +
         qexpr::_CNOT(q0, q1);
}

PROTECT QExpr rotationQMUX(datalist::DataList p, qlist::QList qs, qbit &q,
                           std::vector<double> &thetas) {

  // For the recursive case, the vector thetas (of size 2^N)
  // should be broken up into newthetas0 and newthetas1,
  // both of size 2^(N-1). First, split thetas in half.
  std::vector<double> thetas0;
  std::vector<double> thetas1;
  for (size_t i = 0; i < thetas.size(); i += 2) {
    thetas0.push_back(thetas[i]);
    thetas1.push_back(thetas[i + 1]);
  }
  std::vector<double> newthetas0; // = 0.5 thetas0 + 0.5 thetas1
  std::vector<double> newthetas1; // = 0.5 thetas0 - 0.5 thetas1
  for (int i = 0; i < thetas0.size(); i++) {
    newthetas0.push_back(0.5 * (thetas0[i] + thetas1[i]));
    newthetas1.push_back(0.5 * (thetas0[i] - thetas1[i]));
  }

  return qexpr::cIf(
      qs.size() == 1,
      // base case
      // singleQubitRotation(p, q, thetas[0]),
      rotationQMUX(p, qs[0], q, thetas[0], thetas[1]),

      // recursive case: qs = qs[0] + (qs+1)
      // equivalent to:
      //    qexpr::qIf(qs[0],
      //        rotationQMUX(p, qs+1, q, thetas1),
      //        rotationQMUX(p, qs+1, q, thetas0)
      //    )
      rotationQMUX(p, qs + 1, q, newthetas0) + qexpr::_CNOT(qs[0], q) +
          rotationQMUX(p, qs + 1, q, newthetas1) + qexpr::_CNOT(qs[0], q));
}

// ============================
// N-qubit state preparation //
// ============================

// Returns a QExpr U such that
//      U|psi> = |next_psi> \otimes |0>
// where |next_psi> is the n-qubit state with cth entry r_c e^{i t_c}
// as specified in the proof of Theorem 9 [ref 1].
//
// Input:
//      qs  - a list of qubits of size N
//      q   - the qubit to disentangle
//      psi - a state vector of size 2^{N+1}
QExpr disentangleLastQubit(qlist::QList qs, qbit &q,
                           std::vector<double> &thetas,
                           std::vector<double> &phis) {
  // Use psi to compute two vectors:
  //      thetas - a vector of 2^N rotation parameters with entries theta_c
  //      phis   - a vector of 2^N rotation parameters with entries phi_c

  return rotationQMUX("Y", qs, q, thetas) * rotationQMUX("Z", qs, q, phis);
}

// Returns a QExpr U such that U |psi> = |0..0>
// Input:
//    qs - QList of size N > 0
//    psi - state vector of size 2^N
QExpr disentangleState(qlist::QList qs, std::vector<C> psi, bool debug) {
  // for the recursive case, compute next_psi, thetas, phis
  std::vector<C> next_psi;
  std::vector<double> thetas;
  std::vector<double> phis;

  for (auto it = psi.begin(); it != psi.end();
       std::ranges::advance(it, 2, psi.end())) {
    auto [r_c, t_c, theta_c, phi_c] = braketToSpherical(*it, *(it + 1), debug);

    C psi_entry = r_c * std::exp(-i * t_c / 2.0);
    next_psi.push_back(psi_entry);
    thetas.push_back(theta_c);
    phis.push_back(phi_c);
  }

  return qexpr::cIf(qs.size() == 1,
                    // base case, |qs|=1, |psi|=2
                    disentangle1(qs[0], psi[0], psi[1], debug),

                    // recursive case, qs = qs(1,qs.size()) + qs[0]
                    disentangleLastQubit(qs + 1, qs[0], thetas, phis) +
                        disentangleState(qs + 1, next_psi, debug));
}

QExpr prepareState(qlist::QList qs, std::vector<C> psi, bool debug) {
  return qexpr::map(qexpr::_PrepZ, qs) +
         qexpr::invert(disentangleState(qs, psi, debug));
}

#endif
