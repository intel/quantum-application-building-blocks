//===--- unitary_ops.h ----------------------------*- C++ -*---------------===//
//
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
/// \file unitary_ops.h
/// \brief A library of unitary operations.
///
/// A set of \c quantum_kernel and \c QExpr functions that apply a
/// unitary transformation to a qubit or a unitary transformation to a qubit
/// conditioned by one or more qubits, i.e. U, controlled-U, and cU gates.
///
//===----------------------------------------------------------------------===//

#ifndef UNITARY_OPS_H
#define UNITARY_OPS_H

// Intel(R) Quantum SDK header files
#include <clang/Quantum/quintrinsics.h>
#include <clang/Quantum/qexpr.h>
#include <qexpr_utils.h>

#include <umatrix.h>

/// \var C
/// @brief A type definition for \c std::complex<double>.
using C = std::complex<double>;


//////////////////
// Declarations //
//////////////////

/** 
 * @defgroup UnitaryZYZ Unitary Gates Z-Y-Z
 * These functions apply a unitary gate in terms of rotation gates RZ-RY-RZ
 * according to
 * Barenco, Adriano, et al. "Elementary gates for quantum computation."
 * Physical Review A 52.5 (1995): 3457.
 * or Theorem 4.1 of
 * M. A. Nielsen and I. L. Chuang, Quantum Computation and Quantum Information:
 * 10th Anniversary Edition. Cambridge: Cambridge University Press, 2010.
 *
 * Theorem 4.1 states:
 * Suppose U is a unitary operation on a single qubit, then there exists real
 * numbers
 *    \f$    \alpha, \beta, \gamma, \delta\    \f$
 * such that
 *    \f$    U = \exp^{i\alpha}  RZ(\beta)  RY(\gamma)  RZ(\delta)    \f$
 * The phase factor
 *    \f$    \exp\left\{i\alpha\right\}    \f$
 * is not a physical observable and doesn't affect the probabilities of a single
 * qubit state.
 * @{
 */

/// @brief Apply a unitary gate.
///
/// @param beta the angle of the final RZ rotation applied.
/// @param gamma the angle of the RY rotation applied.
/// @param delta the angle of the first RZ rotation applied.
/// @param target the qbit that the operation is applied to.
quantum_kernel void Uzyz(const double &beta, const double &gamma,
                         const double &delta, qbit &target);

/// @brief Apply a unitary gate.
///
/// @param beta   the angle of the final RZ rotation applied.
/// @param gamma  the angle of the RY rotation applied.
/// @param delta  the angle of the first RZ rotation applied.
/// @param target the qbit that the operation is applied to.
QExpr _Uzyz(qbit &target, const double &beta,
             const double &gamma, const double &delta);

/// @brief Return a QExpr implementing the unitary U on the qubit q.
///
/// @param U A single-qubit (2x2) unitary.
/// @param q The qubit on which to apply \c U.
/// @param debug Perform checks and print additional output for debugging.
QExpr fromUnitary1(const UMatrix<1> &U, qbit &q, bool debug = false);


/** @} */ // end of UnitaryZYZ


/** 
 * @defgroup cUnitaryZYZ Controlled Unitary Gates Z-Y-Z
 * These functions apply a controlled-unitary gate in terms of rotation gates
 * RZ-RY-RZ conditioned on one or more other qubits according to
 * Barenco, Adriano, et al. "Elementary gates for quantum computation."
 * Physical Review A 52.5 (1995): 3457.
 * or Corollary 4.2 of
 * M. A. Nielsen and I. L. Chuang, Quantum Computation and Quantum Information:
 * 10th Anniversary Edition. Cambridge: Cambridge University Press, 2010.
 * Corollary 4.2 states:
 * Suppose
 *    \f$    U    \f$
 * is a unitary gate on a single qubit. Then there exists unitary operators
 *    \f$    A, B, C    \f$
 * on a single qubit such that
 *    \f$    ABC = I    \f$
 * and
 *    \f$    U = \exp\left\{i \alpha \right\} A X B X C    \f$
 * where
 *    \f$    \alpha    \f$
 * is some overall phase factor, and
 *    \f$    X    \f$
 * is the Pauli operator.
 * @{
 */

/// @brief Apply a controlled unitary gate.
///
/// @param ctrl the qbit controlling the operation.
/// @param target the qbit that the operation is applied to.
/// @param beta the angle of the final RZ rotation applied.
/// @param gamma the angle of the RY rotation applied.
/// @param delta the angle of the first RZ rotation applied.
quantum_kernel void cUzyz(qbit &ctrl, qbit &target,
                          const double &beta, const double &gamma,
                          const double &delta);

/// @brief Apply a unitary controlled by two qubits and using one ancilla
///
/// Apply a controlled-controlled unitary gate according to the so-called
/// Toffoli ladder construction.
/// The gate U will be applied if and only if qbits ctrl1 and ctrl2 are set.
///
/// @param ctrl1 qbit variable whose state will determine the first control.
/// @param ctrl2 qbit variable whose state will determine the second control.
/// @param work1 the ancilla qubit that facilitates the Toffoli ladder
/// @param target the qbit variable to apply the gate to.
/// @param beta  the angle of the final RZ rotation applied.
/// @param gamma the angle of the RY rotation applied.
/// @param delta the angle of the first RZ rotation applied.
quantum_kernel void twoControl_U(qbit &ctrl1, qbit &ctrl2, qbit &work1,
                                 qbit &target, const double &beta,
                                 const double &gamma,
                                 const double &delta);

/// @brief Apply a unitary controlled by three qubits and using two ancilla.
///
/// Apply a c-c-c-unitary gate according to the so-called Toffoli ladder
/// construction.
/// The gate U will be applied if and only if qbits ctrl1, ctrl2, and ctrl3 are
/// set.
///
/// @param ctrl1  qbit variable whose state will determine the first control.
/// @param ctrl2  qbit variable whose state will determine the second control.
/// @param ctrl3  qbit variable whose state will determine the third control.
/// @param work1  ancilla qubit that facilitates the Toffoli ladder
/// @param work2  ancilla qubit that facilitates the Toffoli ladder
/// @param target the qbit variable to apply the gate to.
/// @param beta   the angle of the final RZ rotation applied.
/// @param gamma  the angle of the RY rotation applied.
/// @param delta  the angle of the first RZ rotation applied.
quantum_kernel void threeControl_U(qbit &ctrl1, qbit &ctrl2, qbit &ctrl3,
                                   qbit &work1, qbit &work2, qbit &target,
                                   const double &beta, const double &gamma,
                                   const double &delta);

/// @brief Return QExpr of a unitary gate with one or more controls.
///
/// Return a _Uzyz unitary gate with one or more control qubits using
/// Functional Language Extension for Quantum (FLEQ) syntax.
///
/// Here, the _Uzyz function is further transformed by adding a number of
/// controls. The number and identity of the control qubits is determined by the
/// contents of the QList container ctrls.
///
/// @param ctrls  a QList container, each qbit variable is a control.
/// @param target the qbit variable to apply _Uzyz to.
/// @param beta   the angle of the final RZ rotation applied.
/// @param gamma  the angle of the RY rotation applied.
/// @param delta  the angle of the first RZ rotation applied.
QExpr multi_ctrl_Uzyz(qlist::QList ctrls, qbit &target,
                      const double &beta, const double &gamma,
                      const double &delta);


/** @} */ // end of cUnitaryZYZ


////////////////////
// Implementation //
////////////////////

////////////////////////////////////////////
// Single-qubit unitary ZYZ decomposition //
////////////////////////////////////////////

quantum_kernel void Uzyz(const double &beta, const double &gamma,
                         const double &delta, qbit &target){
  RZ(target, delta);
  RY(target, gamma);
  RZ(target, beta);
}


QExpr _Uzyz(qbit &target, const double &beta,
            const double &gamma, const double &delta) {

  return qexpr::_RZ(target, delta) + qexpr::_RY(target, gamma)
         + qexpr::_RZ(target, beta);
}


quantum_kernel void cUzyz(qbit &ctrl, qbit &target, const double &beta,
                          const double &gamma, const double &delta) {
  // C
  double rotation_c = (delta - beta)/2;
  RZ(target, rotation_c);

  // X
  CNOT(ctrl, target);

  // B
  double rotation_bz = -1.0*(delta + beta)/2;
  RZ(target, rotation_bz);

  double rotation_by = gamma / (-2);
  RY(target, rotation_by);

  // X
  CNOT(ctrl, target);

  // A
  double rotation_ay = gamma / 2;
  RY(target, rotation_ay);

  double rotation_az = beta;
  RZ(target, rotation_az);
}


quantum_kernel void twoControl_U(qbit &ctrl1, qbit &ctrl2, qbit &work1,
                                 qbit &target, const double &beta,
                                 const double &gamma,
                                 const double &delta) {
  Toffoli(ctrl1, ctrl2, work1);
  cUzyz(work1, target, beta, gamma, delta);
  Toffoli(ctrl1, ctrl2, work1);
}


quantum_kernel void threeControl_U(qbit &ctrl1, qbit &ctrl2, qbit &ctrl3,
                                   qbit &work1, qbit &work2, qbit &target,
                                   const double &beta, const double &gamma,
                                   const double &delta) {
  Toffoli(ctrl1, ctrl2, work1);
  twoControl_U(ctrl3, work1, work2, target, beta, gamma, delta);
  Toffoli(ctrl1, ctrl2, work1);
}


QExpr multi_ctrl_Uzyz(qlist::QList ctrls, qbit &target,
                      const double &beta, const double &gamma,
                      const double &delta) {

  return qexpr::control(ctrls, _Uzyz(target, beta, gamma, delta));
}

/// @brief Return a tuple of the angles
///       \f$U = \exp\left\{i\phi\right\} RZ(\alpha) RY(\beta) RZ(\gamma)\f$.
///
/// Return a tuple of angles such that the unitary gate defined by the
/// input (2x2) matrix U is equivalent to
/// \f$ U = e^{i \phi} \left[ \begin{matrix} e^{-i\gamma/2 -i\alpha/2} \cos(\beta/2) & -e^{-i\gamma/2 +i\alpha/2} \sin(\beta/2) \\ e^{ i\gamma/2 -i\alpha/2} \sin(\beta/2) &  e^{ i\gamma/2 +i\alpha/2} \cos(\beta/2) \\ \end{matrix} \right] \f$
/// Each row of the matrix is an entry in the outer vector. The entries of the
/// row are the elements of the inner vector.
///
/// @param U A single-qubit (2x2) unitary matrix given as a vector of vectors.
/// @param debug Perform checks and print additional output for debugging.
// Return a tuple of the angles such that
// U = e^{i * phi} RZ(alpha) RY(beta) RZ(gamma).
///
// Return a tuple of angles such that the unitary gate defined by the
// input (2x2) matrix U is equivalent to
// U = e^{i*phi} *
//|e^{-i*gamma/2-i*alpha/2}*cos(beta/2) -e^{-i*gamma/2+i*alpha/2}*sin(beta/2)|
//|e^{ i*gamma/2-i*alpha/2}*sin(\beta/2) e^{i*gamma/2+i*alpha/2}*cos(beta/2) |
std::tuple<double, double, double, double>
unitary1ToZYZDecomposition(const UMatrix<1> &U, bool debug = false) {
  const double EPSILON = 10e-4; // used for determining == 0

  // compute Phi and rescale the remaining elements
  // determinant = e^(i*2*Phi)
  C determinant = U.determinant();
  double Phi = std::atan2(determinant.imag(), determinant.real()) / 2;

  C phase(std::cos(Phi), -std::sin(Phi));
  auto v = phase * U;

  double beta = std::acos(2 * std::norm(v(0, 0)) - 1); // assume |v00| != 0
  if (isEqual(v(0, 0), 0.0, EPSILON)) {
    beta = std::asin(1 - 2 * std::norm(v(0, 1)));
  }

  double alpha_plus_gamma = 2 * std::atan2(v(1, 1).imag() / std::cos(beta / 2),
                                           v(1, 1).real() / std::cos(beta / 2));
  if (isEqual(std::cos(beta / 2), 0.0, EPSILON)) {
    alpha_plus_gamma = 0;
  }

  double alpha_minus_gamma =
      2 * std::atan2(v(1, 0).imag() / std::sin(beta / 2),
                     v(1, 0).real() / std::sin(beta / 2));
  if (isEqual(std::sin(beta / 2), 0.0, EPSILON)) {
    alpha_minus_gamma = 0.0;
  }

  double alpha = (alpha_plus_gamma + alpha_minus_gamma) / 2;
  double gamma = (alpha_plus_gamma - alpha_minus_gamma) / 2;

  if (debug) {
    // Validate that U = e^{i Phi} RZ(alpha) RY(beta) RZ(gamma)
    // C phase(std::cos(Phi), -std::sin(Phi));
    UMatrix<1> RHS =
        std::conj(phase) * RZMatrix(alpha) * RYMatrix(beta) * RZMatrix(gamma);

    if (!isEqual<1>(U, RHS, EPSILON)) {
      std::cout << "== ERROR in unitary1ToZYZDecomposition ==\n";
      std::cout << "U: \n" << U << "\n";
      std::cout << "Phi: " << Phi << "\n";
      std::cout << "alpha: " << alpha << "\n";
      std::cout << "beta: " << beta << "\n";
      std::cout << "gamma: " << gamma << "\n";
      std::cout << "RHS: \n" << RHS << "\n";
      assert(false);
    }
  }

  return {Phi, alpha, beta, gamma};
}


QExpr fromUnitary1(const UMatrix<1> &U, qbit &q, bool debug) {
    auto [Phi, alpha, beta, gamma] = unitary1ToZYZDecomposition(U);
    // Ignore the phase (first tup) because we only work up to phase
    return qexpr::_RZ(q, alpha) + qexpr::_RY(q, beta) + qexpr::_RZ(q, gamma);
}

#endif
