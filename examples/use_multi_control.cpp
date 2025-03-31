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
/// \file use_multi_control.cpp
/// \brief Demonstrate adding controls to gates using each type of
/// quantum function, ``quantum_kernel`` and ``QExpr``.
///
//===----------------------------------------------------------------------===//
// Intel(R) Quantum SDK header files
#include <clang/Quantum/quintrinsics.h>
#include <qexpr_utils.h>
#include <quantum_full_state_simulator_backend.h>

// qabbl/include/
#include <unitary_ops.h>

const int total_qubits = 4;
qbit qubit_register[total_qubits];

const double test_angle = M_PI_2;

/*
Circuit to test Uzyz(), performed on 2 qubits to compare with controlled U case
*/
quantum_kernel void test_U() {
  double beta = test_angle;
  double gamma = test_angle;
  double delta = test_angle;
  for (int i = 0; i < total_qubits; i++) {
    PrepZ(qubit_register[i]);
  }

  X(qubit_register[1]);

  Uzyz(beta, gamma, delta, qubit_register[0]);
}

/*
Circuit to test cUzyz()
*/
quantum_kernel void test_cU() {
  double beta = test_angle;
  double gamma = test_angle;
  double delta = test_angle;
  for (int i = 0; i < total_qubits; i++) {
    PrepZ(qubit_register[i]);
  }
  X(qubit_register[1]);
  cUzyz(qubit_register[1], qubit_register[0], beta, gamma, delta);
}

/*
Circuit to test twoControl_U()
*/
quantum_kernel void test_twoControl_U() {
  double beta = test_angle;
  double gamma = test_angle;
  double delta = test_angle;
  for (int i = 0; i < total_qubits; i++) {
    PrepZ(qubit_register[i]);
  }

  X(qubit_register[1]);
  X(qubit_register[2]);
  twoControl_U(qubit_register[1], qubit_register[2], qubit_register[3],
               qubit_register[0], beta, gamma, delta);
}

/*
Circuit to test _Uzyz().
*/
PROTECT QExpr test_U_fleq() {
  double beta = test_angle;
  double gamma = test_angle;
  double delta = test_angle;

  for (int i = 0; i < total_qubits; i++) {
    PrepZ(qubit_register[i]);
  }

  X(qubit_register[1]);

  return this_as_expr + _Uzyz(qubit_register[0], beta, gamma, delta);
}

/*
Circuit to test multi_ctrl_Uzyz() with a single control.
*/
PROTECT QExpr test_ctrl_U_fleq() {
  double beta = test_angle;
  double gamma = test_angle;
  double delta = test_angle;

  for (int i = 0; i < total_qubits; i++) {
    PrepZ(qubit_register[i]);
  }

  X(qubit_register[1]);

  qlist::QList controls(qubit_register[1]);

  return this_as_expr +
         multi_ctrl_Uzyz(controls, qubit_register[0], beta, gamma, delta);
}

/*
Circuit to test multi_ctrl_Uzyz() with a two controls.
*/
PROTECT QExpr test_ccU_fleq() {
  double beta = test_angle;
  double gamma = test_angle;
  double delta = test_angle;

  for (int i = 0; i < total_qubits; i++) {
    PrepZ(qubit_register[i]);
  }

  X(qubit_register[1]);
  X(qubit_register[2]);

  qlist::QList control1(qubit_register[1]);
  qlist::QList control2(qubit_register[2]);
  qlist::QList controls = control1 + control2;

  return this_as_expr +
         multi_ctrl_Uzyz(controls, qubit_register[0], beta, gamma, delta);
}

/// @cond
int main() {
  iqsdk::IqsConfig settings(total_qubits, "noiseless");
  iqsdk::FullStateSimulator quantum_8086(settings);
  if (iqsdk::QRT_ERROR_SUCCESS != quantum_8086.ready()) {
    return 1;
  }
  // Create a references to qbits we want data from
  std::vector<std::reference_wrapper<qbit>> qids;
  // Case A: output will show the state of the entire qubit_register
  //         to focus only on the target qubit, comment out this loop and
  //         uncomment Case B below
  for (int id = 0; id < 2; ++id) {
    qids.push_back(std::ref(qubit_register[id]));
  }
  // End Case A.

  // Case B: output will show only the target qubit
  // qids.push_back(std::ref(qubit_register[0]));
  // end Case B.

  // build a wildcard string for input to patternToIndices
  std::string wildcard;
  for (int j = 0; j < 2; ++j) {
    wildcard.push_back('?');
  }
  // create a QssIndex for each possible state
  std::vector<iqsdk::QssIndex> all_states =
      iqsdk::QssIndex::patternToIndices(wildcard);

  std::cout << "Running Uzyz()" << std::endl;
  test_U(); // quantum_kernel functions are called directly
  iqsdk::QssMap<double> probability_map =
      quantum_8086.getProbabilities(qids, all_states);
  quantum_8086.displayProbabilities(probability_map);

  std::cout << "Running _Uzyz()" << std::endl;
  qexpr::eval_hold(test_U_fleq()); // QExpr functions are passed to eval_hold
  probability_map = quantum_8086.getProbabilities(qids, all_states);
  quantum_8086.displayProbabilities(probability_map);

  std::cout << "Running controlled-Uzyz()" << std::endl;
  test_cU();
  probability_map = quantum_8086.getProbabilities(qids, all_states);
  quantum_8086.displayProbabilities(probability_map);

  std::cout << "Running controlled-_Uzyz()" << std::endl;
  qexpr::eval_hold(test_ctrl_U_fleq());
  probability_map = quantum_8086.getProbabilities(qids, all_states);
  quantum_8086.displayProbabilities(probability_map);

  std::cout << "Running controlled-controlled-Uzyz()" << std::endl;
  qids.push_back(std::ref(qubit_register[2]));
  wildcard.push_back('?');
  all_states = iqsdk::QssIndex::patternToIndices(wildcard);
  test_twoControl_U();
  probability_map = quantum_8086.getProbabilities(qids, all_states);
  quantum_8086.displayProbabilities(probability_map);

  std::cout << "Running controlled-controlled-_Uzyz()" << std::endl;
  qexpr::eval_hold(test_ccU_fleq());
  probability_map = quantum_8086.getProbabilities(qids, all_states);
  quantum_8086.displayProbabilities(probability_map);

  return 0;
}
/// @endcond
