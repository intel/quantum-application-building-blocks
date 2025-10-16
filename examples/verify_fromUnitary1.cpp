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
/// \file verify_fromUnitary1.cpp
/// \brief Demonstration and test of fromUnitary1() .
///
/// The example runs \c haar_distribution() in a loop, creating
/// a sequence of unitary 2x2 matrices. Each matrix is converted to a quantum
/// algorithm and applied to two reference quantum states using the
/// Intel Quantum Simulator. Those results are checked against the equivalent
/// operations on the matrix elements to within \c 1e-4.
///
/// The conjugate transpose of the matrix is also converted to quantum algorithm
/// and applied to ensure that the original quantum state is recovered to within
/// \c 1e-7.
//===----------------------------------------------------------------------===//
// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/quintrinsics.h>
#include <quantum_full_state_simulator_backend.h>

// qabbl/include/
#include <unitary_ops.h>

// C++ standard library
#include <iostream>

QExpr verifyCompZero(const UMatrix<1> &U, qbit &q) {
  return qexpr::_PrepZ(q) + fromUnitary1(U, q);
}

QExpr verifyCompOne(const UMatrix<1> &U, qbit &q) {
  return qexpr::_PrepZ(q) + qexpr::_X(q) + fromUnitary1(U, q);
}

QExpr formIdentity(UMatrix<1> &U, qbit &q) {

  return fromUnitary1(U.adjoint(), q);
}

/// @cond
int main() {
  double eps = 1e-4; // for comparing probabilities

  // an engine that produces a sequence of pseudo-random values.
  std::mt19937 engine(1234); // fixed seed for deterministic sequence

  // std::random_device rd; // random seed
  // std::mt19937 engine(rd()); // non-deterministic sequence

  iqsdk::IqsConfig settings;
  iqsdk::FullStateSimulator iqs_device(settings);
  iqsdk::QRT_ERROR_T status = iqs_device.ready();
  assert(status == iqsdk::QRT_ERROR_SUCCESS);

  qbit q0;

  std::vector<std::reference_wrapper<qbit>> qids;
  qids.push_back(q0);
  std::vector<iqsdk::QssIndex> bases = iqsdk::QssIndex::patternToIndices("?");

  for (int tries = 0; tries < 100000; tries++) {
    UMatrix<1> matrix_rep = haar_distribution<1>(engine);

    for (int state = 0; state < 2; state++) {
      /* Compare Prob( U | state > ) against U_state0 and U_state1 */
      if (state == 0) {
        qexpr::eval_hold(verifyCompZero(matrix_rep, q0));
      }
      if (state == 1) {
        qexpr::eval_hold(verifyCompOne(matrix_rep, q0));
      }
      iqsdk::QssMap<double> probability_map =
          iqs_device.getProbabilities(qids, bases);
      double measureZero = probability_map[bases.at(0)] -
        std::real(matrix_rep(state, 0) *
                  std::conj(matrix_rep(state, 0)));
      double measureOne = probability_map[bases.at(1)] -
        std::real(matrix_rep(state, 1) *
                  std::conj(matrix_rep(state, 1)));

      if (measureZero > eps) {
        std::cout << "The difference between IQS probabilities and matrix "
                     "elements is too large!\n";
        std::cout << "state: " << state << std::endl;
        std::cout << "the complex-square of the matrix elements\n";
        std::cout << matrix_rep(state,0) * matrix_rep.conjugate()(state,0)
                  << std::endl;
        std::cout << "the IQS-probabilities\n";
        iqs_device.displayProbabilities(probability_map);
        std::cout << "difference in measuring 0 " << measureZero << std::endl;

        return 1;
      } // end measureZero
      if (measureOne > eps) {
        std::cout << "The difference between IQS probabilities and matrix "
                     "elements is too large!\n";
        std::cout << "state: " << state << std::endl;
        std::cout << "the complex-square of the matrix elements\n";
        std::cout << matrix_rep(state,1) * matrix_rep.conjugate()(state,1)
                  << std::endl;
        std::cout << "the IQS-probabilities\n";
        iqs_device.displayProbabilities(probability_map);
        std::cout << "difference in measuring 1 " << measureOne << std::endl;
        std::cout << matrix_rep * matrix_rep.conjugate() << std::endl;
        return 1;
      } // end measureOne
      /* Apply U-dagger and compare to initial probability */
      qexpr::eval_hold(formIdentity(matrix_rep, q0));
      probability_map = iqs_device.getProbabilities(qids, bases);
      double unity_difference = 1.0 - probability_map[bases.at(state)];
      if (unity_difference > eps) {
        std::cout
            << "The difference between U U^dagger and Identity is too large!\n";
        std::cout << "state: " << state << std::endl;
        std::cout << "the complex-square of the matrix elements\n";
        std::cout << matrix_rep * matrix_rep.adjoint() << std::endl;
        std::cout << "Unity difference: " << unity_difference << std::endl;
        iqs_device.displayProbabilities(probability_map);
        return 1;
      }
    } // loop over state
  }   // loop over tries

  return 0;
}
/// @endcond
