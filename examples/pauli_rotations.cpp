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
// This file demonstrates the use of three functions:
// pauliRotation(), pauliMeas(), and pauliPrep(). These operations
// take as input a Pauli string and output a quantum kernel expression
// implementing the multi-qubit rotation, measurement, or preparation specified.
//
//===----------------------------------------------------------------------===//
// Intel(R) Quantum SDK header files
#include <quantum_full_state_simulator_backend.h>

// qabbl/include/
#include <pauli_ops.h>

// C++ standard library
#include <cassert>

int main() {
    iqsdk::IqsConfig iqs_config;
    iqsdk::FullStateSimulator iqs_device(iqs_config);
    iqsdk::QRT_ERROR_T status = iqs_device.ready();
    assert(status == iqsdk::QRT_ERROR_SUCCESS);

    const int N = 12;
    qbit listable(qs,N);


    qexpr::eval_hold(qexpr::printQuantumLogic(
            pauliRotation("0.03 Z0 X1", qs)
    ));


    qexpr::eval_hold(qexpr::printQuantumLogic(
            // By default theta=1
            pauliRotation("X10 Y0 Z3", qs)
    ));
    qexpr::eval_hold(qexpr::printQuantumLogic(
            // implements a global phase
            pauliRotation("10", qs)
    ));
    qexpr::eval_hold(qexpr::printQuantumLogic(
            // Can also specify the angle as an argument
            pauliRotationBy(0.01, "Z0", qs)
    ));
    qexpr::eval_hold(qexpr::printQuantumLogic(
            pauliPrep("Y1 Z0 X1", qs)
    ));
    bool b;
    qexpr::eval_hold(qexpr::printQuantumLogic(
            pauliMeas("Y1 Z0", qs, b)
    ));



    qexpr::eval_hold(qexpr::printQuantumLogic(
            pauliRotations("{0.01 Y1 Z0; 3 X0; Z0}", qs)
    ));
    qexpr::eval_hold(qexpr::printQuantumLogic(
            pauliRotations("0.01 Y1 Z0; 3 X0; Z0", qs)
    ));

    qexpr::eval_hold(qexpr::printQuantumLogic(
            qexpr::mapDataList("{", "}", pauliPrep, ";",
                        "{X0 Y1; Z0}",
                        qs)
    ));




    return 0;
}
