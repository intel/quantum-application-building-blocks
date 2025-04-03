//===--- qalgorithm.h -----------------------------*- C++ -*---------------===//
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
/// \file qalgorithm.h
/// \brief A library of quantum algorithms.
///
///  A collection of implementations for using and applying quantum algorithms.
///
//===----------------------------------------------------------------------===//

#ifndef QALGORITHM_H
#define QALGORITHM_H

// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/qlist.h>


//////////////////
// Declarations //
//////////////////

/** 
 * @defgroup qalgo Quantum Algorithms
 * Tools to apply quantum algorithms to qubits.
 * @{
 */

/// @brief  Quantum Fourier Transform
///
/// Applies quantum Fourier transform to a set of qubits.
///
/// QFT is the quantum analogue of discrete Fourier transform. QFT can be thought
/// of as a change of basis from the computational basis
///    \f$    \ket{x}    \f$
/// to the Fourier basis
///    \f$    \ket{\sim x}    \f$,
/// where
///    \f[
///           \ket{\sim x} =  \frac{1}{\sqrt{2^n}} *
///                           \sum_{0 <= y < 2^n} e^{2 \pi i x y / 2^n} \ket{y}
///    \f]
///
/// @param reg A list of qubits.
///
/// @warning Assumes "little endian" number encoding.
QExpr qft(qlist::QList reg);


/// @brief Quantum Phase Estimation
///
/// Applies quantum phase estimation to compute the eigenvalue of the unitary
/// given as input.
///
/// QPE is a subroutine which takes in any unitary U, and approximates
/// the eigenvalue into a separate register for a prepared eigenstate.
///
/// More specifically, if we can prepare the state
///    \f$    \ket{\psi}    \f$
/// such that if
///    \f$    U \ket{\psi} = e^{2 \pi i \theta} \ket{\psi}    \f$,
/// where
///    \f$    0 \le \theta < 1    \f$ ,
/// then QPE takes the state
///    \f$    \ket{\psi} \ket{0}    \f$
/// to the state
///    \f$    \ket{\psi}{\theta'} + O(e)    \f$,
/// where
///    \f$    |\theta - \theta'| < 1/2^N    \f$
/// for a phase register of size
///    \f$    N    \f$.
///
/// @param reg A list of qubits used for the phase register.
/// @param U the unitary whose phase is to be estimated
///
/// @warning Assumes "little endian" number encoding.
QExpr qpe(qlist::QList reg, QExpr U);

/** @} */ // end of qalgo


////////////////////
// Implementation //
////////////////////

/// @brief Reverse the order of qubits in the register
PROTECT QExpr reverseRegister(qlist::QList reg){
  unsigned sz = reg.size();

  return qexpr::cIf(sz > 1,
            reverseRegister(reg(1,sz-1)) + qexpr::_SWAP(reg[0], reg[sz-1])
          , qexpr::identity()
        );
}

/// @brief The inner loop of QFT.
PROTECT QExpr qftCPhaseLadder(qbit& q, qlist::QList reg, double angle){
  int sz = reg.size();
  return qexpr::cIfTrue(sz > 0,
            qexpr::control(reg[0], qexpr::_RZ(q, angle) + qexpr::global_phase(-angle))
            + qftCPhaseLadder(q, reg + 1, angle / 2.)
  );
}


/// @brief The outer loop of the QFT.
PROTECT QExpr qftHelper(qlist::QList reg){
  int sz = reg.size();
  return qexpr::cIfTrue(sz > 0,
                    qexpr::_H(reg[0])
                    + qftCPhaseLadder(reg[0], reg + 1, M_PI / 2.)
                    + qftHelper(reg + 1)
  );
}

//NOTE: '*' is used here because of the little endian convention
QExpr qft(qlist::QList reg){
  return qftHelper(reg) * reverseRegister(reg);
}


/*
The primary driver of QPE is the control of unitaries U(m)
for 0 <= m < size(reg). These unitaries are just the passed 
unitary, U, to an integer power which is itself a power of 2.
i.e. U(m) = U^(2^m). This allows for a nice recursive form,
U(1) = U,
U(m + 1) = U(m)^2
which we can translate to FLEQ recursion through the following
helper function.
Note, as this is the little endian convention, the least
significant (smallest power) corresponds to the beginning of the
register.
*/
QExpr qpeHelper(qlist::QList reg, QExpr Um){
  unsigned sz = reg.size();
  return qexpr::cIfTrue(sz > 0,
              //control the passed unitary (i.e. this is the mth step)
              //on the first qubit in reg
              qexpr::control(reg[0], Um)
              //recurse by passing reg slicing at the begining and Um squared
              + qpeHelper(reg >> 1, Um^2)
            );
}

// QPE typically start in the all |0> state so we will enforce that in the QExpr
// function:
QExpr qpe(qlist::QList reg, QExpr U){
  return qexpr::map(qexpr::_PrepZ, reg) //prep to all |0> state
         + qexpr::map(qexpr::_H, reg) // create even superposition
         + qpeHelper(reg, U) // apply controlled unitaries
         + ~qft(reg); // apply inverse QFT
}

// We may also want a pure unitary version, where we need to define how it acts,
// on all basis state. For this, we might define a "cyclic" version:
// |psi>|a> -(QPE)-> |psi>|(theta + a)' mod(2^N)> + O(e).

//This simply requires the Preps and Hadamards be replaced with QFT:
QExpr qpe_cyclic(qlist::QList reg, QExpr U){
  return qft(reg) //apply QFT
         + qpeHelper(reg, U) // apply controlled unitaries
         + ~qft(reg); // apply inverse QFT
}

#endif
