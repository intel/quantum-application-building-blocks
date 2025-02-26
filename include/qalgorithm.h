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

#ifndef QALGORITHM_H
#define QALGORITHM_H

// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/qlist.h>

//////////////////
// Declarations //
//////////////////

/*******************************************************************************
 * @brief  Quantum Fourier Transform
 *
 * Return a QExpr that applies quantum fourier transform to a set of qubits.
 * Assumes "big endian" number encoding.
 *
 * @param reg  A list of qubits
 ******************************************************************************/
QExpr qft(qlist::QList reg);

/*******************************************************************************
 * @brief Quantum Phase Estimation
 *
 * Use quantum phase estimation to compute the eigenvalue of the unitary given
 * as input.
 * Assumes "big endian" number encoding.
 *
 * @param reg A list of qubits used for the phase register.
 * @param U the unitary whose phase is to be estimated
 ******************************************************************************/
QExpr qpe(qlist::QList reg, QExpr U);


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

//NOTE: '+' is used here because of the big endian convention
QExpr qft(qlist::QList reg){
  return qftHelper(reg) + reverseRegister(reg);
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
Note, as this is the big endian convention, the most 
significant (largest power) corresponds to the beginning of the
register.
*/
QExpr qpeHelper(qlist::QList reg, QExpr Um){
  unsigned sz = reg.size();
  return qexpr::cIfTrue(sz > 0,
              //control the passed unitary (i.e. this is the mth step)
              //on the first qubit in reg
              qexpr::control(reg[sz - 1], Um)
              //recurse by passing reg slicing at the begining and Um squared
              + qpeHelper(reg << 1, Um^2)
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
