//===--- pauli_ops.h ------------------------------*- C++ -*---------------===//
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
/// \file pauli_ops.h
/// \brief A library of operations for Pauli strings.
///
/// A set of operations that take as input a Pauli string and output a
/// quantum kernel expression implementing the multi-qubit rotation,
/// measurement, or preparation specified.
///
//===----------------------------------------------------------------------===//

#ifndef PAULI_OPS_H
#define PAULI_OPS_H

// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/qlist.h>
#include <clang/Quantum/datalist.h>
#include <qexpr_utils.h>


//////////////////
// Declarations //
//////////////////

/** 
 * @defgroup pauliStrings Pauli strings
 * Functions to use Pauli strings to specify qubit preparation, rotations, and
 * measurements.
 * @{
 */


/// @brief      Returns a QExpr that implements the multi-qubit
///             rotation specified by d
/// @param d    A DataList of the form 'P' or 'theta P' where
///             theta is a double and P is a Pauli string.
///             that matches the regular expression
///                     ((X|Y|Z)(0|1|2|3|4|5|6|7|8|9)+)*
///             Whitespace is also allowed.
///             The integers following each (X|Y|Z) should range
///             from 0 to reg.size()-1.
QExpr pauliRotation(datalist::DataList d, qlist::QList qs);

/** @} */ // end of pauliStrings


/// @brief      Returns a QExpr that implements the multi-qubit
///             rotation specified by d around angle theta
PROTECT QExpr pauliRotationBy(double theta, datalist::DataList d, qlist::QList qs);

/// @brief      Returns a QExpr that implements one or more multi-qubit
///             rotations specified by d
/// @param d    A DataList of the form '{ d1 ; ... ; dn }' or just 'd1 ; ... ; dn' where
///             each d is a valid DataList input to 'pauliRotation'.
QExpr pauliRotations(datalist::DataList d, qlist::QList qs);

/** @addtogroup pauliStrings
 *  
 *  @{
 */

/// @brief    Returns a QExpr that implements a multi-qubit Pauli preparation
///           specified by d.
/// @param d  A Pauli string in the form of a DataList, where all indices
///           range from 0 to reg.size()-1.
PROTECT QExpr pauliPrep(datalist::DataList d, qlist::QList qs);


/// @brief    Returns a QExpr that implements the multi-qubit Pauli measurement
///           specified by d, writing the result to b.
/// @param d  A Pauli string in the form of a DataList, where all indices
///           range from 0 to reg.size()-1.
PROTECT QExpr pauliMeas(datalist::DataList d, qlist::QList qs, bool& b);

/** @} */ // end of addtogroup pauliStrings

/////////////////////////////////
// Multi-Qubit Pauli Rotations //
/////////////////////////////////


/// @brief      Returns a QExpr e such that
///                 Rot(Z(q1)p2(q2))
///             is equivalent to
///                 e + Rot(p2(q2)) + invert(e)
/// @param p2   A QList of length 1 equal to "X", "Y", or "Z"
QExpr reduceSupportZ(qbit& q1, datalist::DataList p2, qbit& q2) {
    return  qexpr::cIf(p2 == datalist::DataList("X"), qexpr::_CZ(q1, q2),
            qexpr::cIf(p2 == datalist::DataList("Y"), qexpr::_CNOT(q1, q2),
            qexpr::cIf(p2 == datalist::DataList("Z"), qexpr::_CNOT(q1, q2),
            qexpr::exitAtCompile(
                "reduceSupportZ: Expected X, Y, or Z; got "
                + p2
                )
    )));
}

/// @brief      Returns a QExpr e such that
///                 Rot(p1(q1)p2(q2))
///             is equivalent to
///                 e + Rot(p2(q2)) + invert(e)
/// @param p1   A QList of length 1 equal to "X", "Y", or "Z"
/// @param p2   A QList of length 1 equal to "X", "Y", or "Z"
QExpr reduceSupport(datalist::DataList p1, qbit& q1, datalist::DataList p2, qbit& q2) {
    return  qexpr::cIf(p1 == datalist::DataList("Z"),
                    reduceSupportZ(q1,p2,q2),
            qexpr::cIf(p1 == datalist::DataList("X"),
                    // H X H = Z
                    qexpr::_H(q1) * reduceSupportZ(q1,p2,q2),
            qexpr::cIf(p1 == datalist::DataList("Y"),
                    // H Sdag Y S H = Z
                    qexpr::_S(q1) * qexpr::_H(q1)
                    * reduceSupportZ(q1, p2, q2),
            qexpr::exitAtCompile(
                "reduceSupport: Expected X, Y, or Z; got "
                + p1
                )
            )));
}

/// @brief        Returns a QExpr that implements a single-qubit
///               rotation around the Pauli axis specified by p
///
/// @param p      A \c DataList of length 1 equal to "X", "Y", or "Z" to specify
///               the rotation axis.
/// @param q      The \c qbit variable the rotation operates on.
/// @param theta  The rotation angle about the specified axis.
QExpr singleQubitRotation(datalist::DataList p, qbit& q, double theta) {
    return  qexpr::cIf(p == datalist::DataList("X"), qexpr::_RX(q,theta),
            qexpr::cIf(p == datalist::DataList("Y"), qexpr::_RY(q,theta),
            qexpr::cIf(p == datalist::DataList("Z"), qexpr::_RZ(q,theta),
            qexpr::exitAtCompile(
                "singleQubitRotation: expected X, Y, or Z; got "
                + p
            ))));
}



/// @brief          Returns a QExpr that implements the multi-qubit
///                 rotation `Rot(P,theta)`, where P is a Pauli string
///                 indexing into the QList qs formed by combining
///                 `p_idx` and the Pauli string `d`
///
/// @param p        A DataList of length 1 equal to "X", "Y", or "Z"
/// @param d        A DataList that matches the regular expression
///                     ((X|Y|Z)(0|1|2|3|4|5|6|7|8|9)+)*
///                 Whitespace is also allowed.
///                 The integers following each Pauli
///                 should be valid indices in qs.
PROTECT QExpr pauliRotationHelper(  datalist::DataList p, int idx,
                            datalist::DataList d, double theta,
                            qlist::QList qs) {

    auto p2     = d(0,1);
    auto dIdx2  = d.next_block("0123456789");
    int  idx2   = dIdx2.to_int();
    auto dRest  = d.after_next(dIdx2).trim();


    return qexpr::cIf(d.empty(),
            singleQubitRotation(p, qs[idx], theta),
            qexpr::conjugate(
                reduceSupport(p, qs[idx], p2, qs[idx2]),
                pauliRotationHelper(p2, idx2, dRest, theta, qs)
                )
            );
}


PROTECT QExpr pauliRotationBy(double theta, datalist::DataList d, qlist::QList qs) {
    auto p      = d(0,1);
    auto dIdx   = d.next_block("0123456789");
    auto dRest  = d.after_next(dIdx).trim();

    return pauliRotationHelper(p, dIdx.to_int(), dRest, theta, qs);
}

PROTECT QExpr pauliRotation(datalist::DataList d, qlist::QList qs) {
    // Trim off any whitespaces
    auto dTrim   = d.trim();
    // Trim the angle from the Pauli string
    auto split   = dTrim.find_any("XYZ");
    double theta = dTrim(0,split).to_double();
    datalist::DataList pauli = dTrim(split, dTrim.size());

    return  qexpr::cIf(split == 0,
                        pauliRotationBy(1, pauli, qs),
            qexpr::cIf(split == d.size(),
                        qexpr::global_phase(theta),
            pauliRotationBy(theta, pauli, qs)
    ));
}

QExpr pauliRotations(datalist::DataList d, qlist::QList qs) {
   return qexpr::mapDataList("{", "}", pauliRotation, ";", d, qs);
}

/////////////////////////////
// Multi-qubit preparation //
/////////////////////////////

QExpr singleQubitPrep(datalist::DataList p, qbit& q) {
    return  qexpr::cIf(p == datalist::DataList("X"), qexpr::_PrepX(q),
            qexpr::cIf(p == datalist::DataList("Y"), qexpr::_PrepY(q),
            qexpr::cIf(p == datalist::DataList("Z"), qexpr::_PrepZ(q),
            qexpr::exitAtCompile(
                "singleQubitPrep: expected X, Y, or Z; got "
                + p
            ))));
}

PROTECT QExpr pauliPrepHelper(datalist::DataList p1, int idx1, datalist::DataList d, qlist::QList qs) {
    auto p2     = d(0,1);
    auto dIdx2  = d.next_block("0123456789");
    int  idx2   = dIdx2.to_int();
    auto dRest  = d.after_next(dIdx2).trim();

    return qexpr::cIf(d.empty(),
                singleQubitPrep(p1, qs[idx1]),
                qexpr::conjugate(
                    reduceSupport(p1, qs[idx1], p2, qs[idx2]),
                    pauliPrepHelper(p2, idx2, dRest, qs)
                )
    );
}

PROTECT QExpr pauliPrep(datalist::DataList d, qlist::QList qs) {
    auto p     = d(0,1);
    auto dIdx  = d.next_block("0123456789");
    auto dRest = d.after_next(dIdx).trim();

    return qexpr::cIf(d.empty(),
              qexpr::exitAtCompile("pauliPrep: Expected a non-empty DataList"),
              pauliPrepHelper(p, dIdx.to_int(), dRest, qs)
    );
}


/////////////////////////////
// Multi-qubit Measurement //
/////////////////////////////

QExpr singleQubitMeas(datalist::DataList p, qbit& q, bool &b) {
    return  qexpr::cIf(p == datalist::DataList("X"), qexpr::_MeasX(q,b),
            qexpr::cIf(p == datalist::DataList("Y"), qexpr::_MeasY(q,b),
            qexpr::cIf(p == datalist::DataList("Z"), qexpr::_MeasZ(q,b),
            qexpr::exitAtCompile(
                "singleQubitMeas: expected X, Y, or Z; got "
                + p
            ))));
}



PROTECT QExpr pauliMeas(datalist::DataList p1, int idx1, datalist::DataList d, qlist::QList qs, bool &b) {
    auto p2    = d(0,1);
    auto dIdx2 = d.next_block("0123456789");
    int  idx2  = dIdx2.to_int();
    auto dRest = d.after_next(dIdx2).trim();

    return qexpr::cIf(d.empty(),
                singleQubitMeas(p1, qs[idx1], b),
                qexpr::conjugate(
                    reduceSupport(p1, qs[idx1], p2, qs[idx2]),
                    pauliMeas(p2, idx2, dRest, qs, b)
                )
    );
}


PROTECT QExpr pauliMeas(datalist::DataList d, qlist::QList qs, bool &b) {
    auto p     = d(0,1);
    auto dIdx  = d.next_block("0123456789");
    auto dRest = d.after_next(dIdx).trim();

    return qexpr::cIf(d.empty(),
              qexpr::exitAtCompile("pauliMeas: Expected a non-empty DataList"),
              pauliMeas(p, dIdx.to_int(), dRest, qs, b)
    );
}

#endif
