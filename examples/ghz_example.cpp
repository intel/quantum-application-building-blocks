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
// Implements the Greenberger–Horne–Zeilinger (GHZ) State for an arbitrary
// number of qubits, preparing the state |0...0>+|1...1>
//
//===----------------------------------------------------------------------===//
// Intel(R) Quantum SDK header files
#include <clang/Quantum/quintrinsics.h>
#include <quantum_full_state_simulator_backend.h>
#include <clang/Quantum/qexpr.h>

// qabbl/include/
#include <qstates.h>

// C++ standard library
#include <iostream>
#include <cassert>
#include <vector>

int main() {
    iqsdk::IqsConfig iqs_config;
    iqsdk::FullStateSimulator iqs_device(iqs_config);
    iqsdk::QRT_ERROR_T status = iqs_device.ready();
    assert(status == iqsdk::QRT_ERROR_SUCCESS);

    const int N = 10;
    qbit listable(q, N);

    qexpr::eval_hold(ghz(q));

    // Print out amplitudes
    std::cout << "------- " << N << " qubit GHZ state -------" << std::endl;
    iqsdk::QssIndex zero_vector (std::string(N, '0'));
    iqsdk::QssIndex one_vector  (std::string(N, '1'));

    auto qbit_refs = to_ref_wrappers(q);
    auto amplitude_map = iqs_device.getAmplitudes(qbit_refs,
                                                  {zero_vector, one_vector});
    iqsdk::FullStateSimulator::displayAmplitudes(amplitude_map);

    return 0;
}
