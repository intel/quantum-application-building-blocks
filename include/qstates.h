//===--- qstates.h --------------------------------*- C++ -*---------------===//
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
/// \file qstates.h
/// \brief A library of quantum states.
///
///  A collection of implementations for setting qubits into specific
///  quantum states.
///
//===----------------------------------------------------------------------===//

#ifndef QSTATES_H
#define QSTATES_H

// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <qexpr_utils.h>

//////////////////
// Declarations //
//////////////////

/** 
 * @addtogroup qstates Quantum States
 * Functions representing quantum states with specific identities and
 * extrapolations.
 * @{
 */

/// @brief Implements the Greenberger–Horne–Zeilinger (GHZ) State
///
/// Prepares a Greenberger–Horne–Zeilinger (GHZ) State for on the input qubits
/// for any number of qubits, i.e. preparing the state
///    \f$    \frac{1}{\sqrt{2}} \{\ket{0 \ldots 0} + \ket{1 \dots 1} \}    \f$
///
/// @param qs A list of qubits to create the GHZ state on
PROTECT QExpr ghz(qlist::QList qs);

/** @} */ // end of qstates

////////////////////
// Implementation //
////////////////////

/////////////////////////////////////////////
// Greenberger–Horne–Zeilinger (GHZ) State //
/////////////////////////////////////////////

PROTECT QExpr ghz(qlist::QList qs) {
    int len = qs.size();
    return (
        // Initialize qs to |0..0>.
        qexpr::map(qexpr::_PrepZ,qs)
        // Prepare the first qubit in state |+>.
        + qexpr::_H(qs[0])
        // Entangle qs[i] with qs[i+1] by mapping CNOT over the QList.
        + qexpr::map(qexpr::_CNOT, qs(0,len-1), qs(1,len)));
        // As an example to understand how this invocation of map works,
        // consider the case where qs has length 3. Then
        //    qs(0,len-1) = { qs[0], qs[1] }
        //    qs(1,len)   = { qs[1], qs[2] }
        // Then this clause will apply _CNOT to each column vector:
        //    _CNOT(qs[0], qs[1]) + _CNOT(qs[1], qs[2])
}

#endif
