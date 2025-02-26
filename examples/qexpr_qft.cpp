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
//
// Demonstrates the quantum fourier transform (QFT).
//
// QFT can be thought of as a change of basis from the computational basis {|x>}
// to the Fourier basis {|~x>}, where
//
// |~x> = 1/sqrt(2^n) sum_{0 <= y < 2^n}( e^{2 pi i x y / 2^n} |y>)
//
// In this file we test QFT by preparing a state |~x> in the Fourier basis
// and applying the inverse QFT to it, which should produce |x> in the computational
// basis.
//
//===----------------------------------------------------------------------===//
// Intel(R) Quantum SDK header files
#include <clang/Quantum/quintrinsics.h>
#include <clang/Quantum/qexpr.h>
#include <quantum_full_state_simulator_backend.h>
#include <qexpr_utils.h>

// qabbl/include/
#include <qalgorithm.h>

// C++ standard library
#include <iostream>
#include <cassert>
#include <vector>


///////////////
/// Testing ///
///////////////

/// @brief  Single-qubit phase gate that implements the unitary
///              (  1   0           )
///              (  0   e^{i theta} )
QExpr phaseGate(qbit& q, double theta) {
    return qexpr::global_phase(-theta) * qexpr::_RZ(q,theta);
}

/// @brief Component to construct the n-qubit fourier basis element |~x> .
///        Applies the appropriate phase gate at qubit q if it appears
///        at index idx.
///
/// @param q    A qubit occurring at index idx in the length-n qubit register
/// @param idx  A qubit index 0 <= idx < n
/// @param x    The target fourier basis element, satisfying 0 <= x < 2^{n-1}
/// @param n    The number of qubits in the total array
PROTECT QExpr fourierPhaseGateAt(qbit& q, int idx, int x, int n) {
    double theta = x * M_PI / pow(2, idx);
    return phaseGate(q, theta);
}


/// @brief      Prepare the n-qubit fourier basis element |~x>.
/// @param qs   An array of qubits of length n
/// @param x    An integer 0 <= x < 2^{n-1}
/// @return
QExpr fourierBasis(qlist::QList qs, int x) {
    return qexpr::map(qexpr::_PrepZ, qs)
            + qexpr::map(qexpr::_H, qs)
            + qexpr::mapWithIndex(fourierPhaseGateAt, qs, x, qs.size());
}


int main() {
    iqsdk::IqsConfig iqs_config;
    iqsdk::FullStateSimulator iqs_device(iqs_config);
    iqsdk::QRT_ERROR_T status = iqs_device.ready();
    assert(status == iqsdk::QRT_ERROR_SUCCESS);


    const int N = 4;
    qbit listable(qs, N);


    // QssIndex's are convertible to integers representing their basis elements.
    iqsdk::QssIndex compBasisIndex("|1011>");
    std::cout << "\n" << "Preparing a fourier basis state corresponding to "
              << compBasisIndex << "\n";
    qexpr::eval_hold(fourierBasis(qs, compBasisIndex) // prepare Fourier basis
                    + -qft(qs) // apply inverse QFT
                    );

    std::cout << "After applying inverse QFT, expect the computational basis element "
              << compBasisIndex << "\n";
    auto qbit_refs = to_ref_wrappers(qs);
    auto probs = iqs_device.getProbabilities(qbit_refs, {}, 0.1);
    iqsdk::FullStateSimulator::displayProbabilities(probs);

    return 0;
}
