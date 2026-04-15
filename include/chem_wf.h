//===--- chem_wf.h ----------------------------------------------*- C++ -*-===//
//
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
/// \file chem_wf.h
/// \brief A library of QExpr functions that implement chemistry wave functions.
///
//===----------------------------------------------------------------------===//

#ifndef CHEM_WF_H
#define CHEM_WF_H

/// @brief Generically prepare all the qubits in the input QList to the zero
///        computational state.
QExpr prep_z_all(qlist::QList qubit_reg) {
    return qexpr::map(qexpr::_PrepZ, qubit_reg);
}

/** 
 * @defgroup HEAWF Hardware Efficient Ansatz
 * These QExpr functions are combined to create the quantum circuit for the
 * Hardware Efficient Ansatz published in
 * Abhinav Kandala, Antonio Mezzacapo, Kristan Temme, Maika Takita,
 * Markus Brink, Jerry M. Chow, and Jay M. Gambetta,
 * Hardware-efficient variational quantum eigensolver for small molecules and
 * quantum magnets. Nature 549, 242 (2017).
 * @{
 */
/// @brief Set the state of the qubit with Ry-Rz-Ry sequence of rotations.
/// @param q      The qubit to perform the instructions on.
/// @param angles A pointer to an array element that is the first of 3 angles
///               needed for 3 rotations on the qubit.
QExpr rotate_qubit(qbit& q, double* angles) {
    return qexpr::_RY(q, *angles) + qexpr::_RZ(q, *(angles + 1)) + qexpr::_RY(q, *(angles + 2));
}

/// @brief Recursively apply rotate_qubit to the first qubit in the input QList
/// @param qubit_reg The list of qubits to rotate, starting at the first entry.
/// @param param_reg A pointer to an array element that is the angle of the
///                  first rotation for the first qubit
QExpr write_to_qubit(qlist::QList qubit_reg, double* param_reg) {
    qlist::QList front = qubit_reg(0,1);
    qlist::QList rest = qubit_reg >> 1;

    return qexpr::cIf(qubit_reg.size() == 1, rotate_qubit(qubit_reg[0], param_reg), rotate_qubit(front[0], param_reg) + write_to_qubit(rest, param_reg + 3) );
}

/// @brief Apply a ladder of CNOT gates to the qubits in the input QList.
/// @param qubit_reg The QList of qubits to apply a sequence of CNOT gates to.
QExpr entangle_qubits(qlist::QList qubit_reg) {
    qlist::QList rest = qubit_reg >> 1;

    return qexpr::cIf(qubit_reg.size() == 2, qexpr::_CNOT(qubit_reg[0], qubit_reg[1]), qexpr::_CNOT(qubit_reg[0], qubit_reg[1]) + entangle_qubits(rest) );
}

/// @brief Apply rotations to set the state of each qubit in the input QList and
///        a then apply a ladder of CNOT gates.
/// @param qubit_reg The qubits to apply operations to.
/// @param param_reg A pointer to an array element that begin the sequence of
///        rotations for all the qubits in qubit_reg.
QExpr one_layer(qlist::QList qubit_reg, double* param_reg) {
    return write_to_qubit(qubit_reg, param_reg) + entangle_qubits(qubit_reg);
}

/// @brief Recursively apply a rotations and CNOT gates in layers.
/// @param qubit_reg The qubits to prepare the layers upon.
/// @param_reg       A pointer to the first array element describing the
///                  3 rotation angles for all qubits in all layers.
/// @param layers    The number of repeated sequences of rotation-entanglement
///                  to apply.
QExpr ansatz_layers(qlist::QList qubit_reg, double* param_reg, unsigned layers) {
    return qexpr::cIf(layers == 1, one_layer(qubit_reg, param_reg), one_layer(qubit_reg, param_reg) + ansatz_layers(qubit_reg, param_reg + 3 * qubit_reg.size(), layers - 1) );
}

/// @brief Reset and prepare the qubits in the hardware efficient ansatz.
/// @param qubit_reg The qubits to prepare the layers upon.
/// @param_reg       A pointer to the first array element describing the
///                  3 rotation angles for all qubits in all layers.
/// @param layers    The number of repeated sequences of rotation-entanglement
///                  to apply.
QExpr hardware_efficient_ansatz(qlist::QList qubit_reg, double* param_reg, unsigned layers) {
    return prep_z_all(qubit_reg) + ansatz_layers(qubit_reg, param_reg, layers);
}

/** @} */ // end of HEAWF
// end Quantum expressions for Hardware Efficient Ansatz
/////////////////////////////////

#endif // CHEM_WF_H
