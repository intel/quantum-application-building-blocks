//===----------------------------------------------------------------------===//
//
// Copyright (C) 2025 Intel Corporation.
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
/// \file validate_prepareState.cpp
/// \brief Demonstration and test of prepareState<N>().
///
/// See "state_prep.h" for detailed descriptions.
///
//===----------------------------------------------------------------------===//
// C++ standard library
#include <cassert>
#include <iostream>
#include <iterator>
#include <vector>

// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/quintrinsics.h>
#include <qexpr_utils.h>
#include <qrt_indexing.hpp>
#include <quantum_full_state_simulator_backend.h>

// qabbl/include/
#include <state_prep.h>
#include <umatrix.h>
#include <unitary_ops.h>
#include <validation_tools.h>


using namespace qlist;
using namespace iqsdk;
using namespace qexpr;


/// @brief Test a sub-routine of prepareState.
template <int N>
bool testRotationQMUX(FullStateSimulator &iqs_device,
                      std::mt19937 &generator, int num_tries = 1,
                      bool debug = true) {
  qbit listable(qs, N);
  qbit q;

  std::vector<double> thetas(powerOfTwo(N));
  std::normal_distribution<double> distribution(1.0, 1.0);

  for (int test = 0; test < num_tries; test++) {
    for (int i = 0; i < thetas.size(); i++) {
      thetas[i] = distribution(generator) / sqrt(2);
      // thetas[i] = i + 0.5;
    }

    UMatrix<N + 1> U = RYMUXMatrix<N>(thetas);

    bool b = implementsUMatrix<N + 1>(
        iqs_device, QList(q) + qs, rotationQMUX<N>("Y", qs, q, thetas, debug),
        U, debug);

    if (!b) {
      std::cout << "Failed testRotationQMUX (1) " << N << "\n";
      assert(false);
    }

    // Also check rotationQMUX when q comes after qs, as in CSD
    // Only for RY
    UMatrix<N + 1> RHS = cosSinMatrix<N>(thetas);

    b = implementsUMatrix<N + 1>(iqs_device, qs + QList(q),
                                 rotationQMUX<N>("Y", qs, q, thetas), RHS,
                                 debug);
    if (!b) {
      std::cout << "Failed testRotationQMUX (2) " << N << "\n";
      assert(false);
    }
  }
  std::cout << "Passed testRotationQMUX for N = " << N << " with " << num_tries
            << " test(s).\n";
  return true;
}


/* This function is left as historical record for the time being; but it is
   commented out so it doesn't impact compile times.

   Each QExpr was expanded so that the qubit state could be checked between
   every transformation step.
*/
// bool debugRotationQMUX(FullStateSimulator &iqs_device) {

//   const int N = 1;

//   QssMap<C> phi;
//   phi[QssIndex("|00>")] = -0.1965 + 0.577 * i;
//   phi[QssIndex("|10>")] = -0.07261 + 0.7028 * i;
//   phi[QssIndex("|01>")] = -0.1155 + 0.3393 * i;

//   QssMap<C> psi;
//   psi[QssIndex("|00>")] = -0.1965 + 0.577 * i;
//   psi[QssIndex("|01>")] = -0.07261 + 0.7028 * i;
//   psi[QssIndex("|10>")] = -0.1155 + 0.3393 * i;

//   std::vector<double> thetas{1.191, 1.95};

//   // debugRotationQMUX2(iqs_device, qssMapToEigen<N+1>(psi), thetas);

//   qbit listable(qs, N);
//   qbit q;

//   assert(isEqual<N + 1>(cosSinMatrix<N>(thetas),
//                         SWAPMatrix * RYMUXMatrix<N>(thetas) * SWAPMatrix));
//   std::cout << "Success (1/3)\n";

//   bool b = implementsUMatrixOn<N + 1>(
//       iqs_device, QList(q) + qs, qssMapToEigen<N + 1>(phi),
//       rotationQMUX<N>("Y", qs, q, thetas, true), RYMUXMatrix<N>(thetas), true);

//   if (!b) {
//     std::cout << "Failed debugRotationQMUX (2)\n";
//     assert(false);
//   }
//   std::cout << "Success (2/3)\n";

//   b = implementsUMatrixOn<N + 1>(
//       iqs_device, qs + QList(q), qssMapToEigen<N + 1>(psi),
//       rotationQMUX<N>("Y", qs, q, thetas, true), cosSinMatrix<N>(thetas), true);

//   if (!b) {
//     std::cout << "Failed debugRotationQMUX (3)\n";
//     assert(false);
//   }
//   std::cout << "Success (3/3)\n";
//   return true;
// }


/* This function is left as historical record for the time being; but it is
   commented out so it doesn't impact compile times.


   Each QExpr was expanded so that the qubit state could be checked between
   every transformation step.
*/
// bool debugRotationQMUX2(FullStateSimulator &iqs_device,
//                         StateVector<2> startingState,
//                         std::vector<double> thetas) {
//   const int N = 1;
//   bool debug = true;
//   assert(thetas.size() == powerOfTwo(N));

//   qbit listable(qs, N);
//   qbit q;
//   auto qsRefs = to_ref_wrappers(qs + QList(q));

//   std::cout << "thetas: ";
//   prettyVector(thetas);
//   std::cout << "\n";

//   eval_hold(prepareState<N + 1>(qs + QList(q), startingState));
//   StateVector<N + 1> sim = startingState;
//   assertState<N + 1>(iqs_device, qsRefs, sim, true);

//   std::cout << "Done with prepareState\n";

//   // The overall unitary being implemented should be:
//   //      RYMUXMatrix(thetas) = directSum(RY(theta0), RY(theta1)).
//   // This should be equal to
//   //      CNOT_{q0, q}
//   //      I ⊗ RY((theta0-theta1) / 2)
//   //      CNOT_{q0, q}
//   //      I ⊗ RY((theta0+theta1) / 2)

//   std::cout << "RYMatrix(+):\n"
//             << RYMatrix((thetas[0] + thetas[1]) / 2.0) << "\n";
//   std::cout << "RYMatrix(-):\n"
//             << RYMatrix((thetas[0] - thetas[1]) / 2.0) << "\n";
//   UMatrix<2> kronPlus = kron<1, 1>(RYMatrix((thetas[0] + thetas[1]) / 2.0),
//                                    UMatrix<1>::Identity());
//   std::cout << "kron+:\n" << kronPlus << "\n";
//   UMatrix<2> kronMinus = kron<1, 1>(RYMatrix((thetas[0] - thetas[1]) / 2.0),
//                                     UMatrix<1>::Identity());
//   std::cout << "kron-:\n" << kronMinus << "\n";

//   UMatrix<2> RHS = SWAPMatrix * CNOTMatrix * SWAPMatrix * kronMinus *
//                    SWAPMatrix * CNOTMatrix * SWAPMatrix * kronPlus;

//   assert(isEqual<2>(cosSinMatrix<N>(thetas), RHS));
//   std::cout << "Done with decomposition\n";

//   eval_hold(_RY(q, (thetas[0] + thetas[1]) / 2.0));
//   sim = kronPlus * sim;
//   assertState<N + 1>(iqs_device, qsRefs, sim);

//   eval_hold(_CNOT(qs[N - 1], q));
//   sim = SWAPMatrix * CNOTMatrix * SWAPMatrix * sim;
//   assertState<N + 1>(iqs_device, qsRefs, sim);

//   eval_hold(singleQubitRotation("Y", q, (thetas[0] - thetas[1]) / 2.0));
//   sim = kronMinus * sim;
//   assertState<N + 1>(iqs_device, qsRefs, sim);

//   eval_hold(_CNOT(qs[N - 1], q));
//   sim = SWAPMatrix * CNOTMatrix * SWAPMatrix * sim;
//   assertState<N + 1>(iqs_device, qsRefs, sim);

//   assertState<N + 1>(iqs_device, qsRefs,
//                      cosSinMatrix<N>(thetas) * startingState);
//   return true;
// }


/// @brief Test a sub-routine of prepareState
template <int N>
bool testDisentangleState(FullStateSimulator &iqs_device,
                          bool debug = false) {

  qbit listable(qs, N);
  auto qsRefs = to_ref_wrappers(qs);
  StateVector<N> psi0 = indexToEigen<N>(QssIndex(N, 0));

  for (int j = 0; j < powerOfTwo(N); j++) {
    QssIndex idx(N, j);
    StateVector<N> psi = indexToEigen<N>(idx);

    eval_hold(prepareQssIndex<N>(qs, idx));

    auto qsRefs = to_ref_wrappers(qs);
    auto probs = iqs_device.getProbabilities(qsRefs, {});
    if (debug) {
      std::cout << "Prepared state (for testDisentangleState): ";
      FullStateSimulator::displayProbabilities(probs);
    }

    if (!implementsOutputState<N>(
            iqs_device, qs, disentangleState<N>(qs, psi, debug), psi0, debug)) {
      std::cout << "Failed testDisentangleState on " << idx << "\n";
      assert(false);
    }
  }

  std::cout << "Passed testDisentangleState for N=" << N << ".\n";
  return true;
}


// Check that
// (1) the QExpr implementation of the Shende-Bullock-Markov algorithm to
// prepare a qubit register in any state up to a global phase successfully
// prepares all basis states and num_tries pseudo-random states drawn from the
// Haar measure.
template <int N>
bool testPrepareState(FullStateSimulator &iqs_device,
                      std::mt19937 &generator, int num_tries = 1,
                      bool debug = true) {
  qbit listable(qs, N);

  for (int j = 0; j < powerOfTwo(N); j++) {
    if (debug) {
      std::cout << "\nprepareState for N=" << N << " qubits; basis state test #"
                << j + 1 << std::endl;
    }
    QssIndex idx(N, j);
    StateVector<N> psi = indexToEigen<N>(idx);

    if (!implementsOutputState<N>(
            iqs_device, qs, prepareState<N>(qs, psi, debug), psi, debug)) {
      std::cout << "Failed testPrepareState with basis state psi = " << idx
                << "\n";
      assert(false);
    }
  }

  for (int i = 0; i < num_tries; i++) {
    if (debug) {
      std::cout << "\nprepareState for N=" << N
                << " qubits; haar random state test #" << i + 1 << std::endl;
    }

    // generate output state
    StateVector<N> psi = haar_state_distribution<N>(generator);

    if (!implementsOutputState<N>(
            iqs_device, qs, prepareState<N>(qs, psi, debug), psi, debug)) {
      std::cout << "Failed testPrepareState with random state\n";
      assert(false);
    }
  }

  std::cout << "Passed testPrepareState for N = " << N << " with " << num_tries
            << " test(s).\n";
  return true;
}


/* This function is left as historical record for the time being; but it is
   commented out so it doesn't impact compile times.


   Each QExpr was expanded so that the qubit state could be checked between
   every transformation step.
*/
// template <int N>
// bool debugPrepareState(FullStateSimulator &iqs_device,
//                        std::mt19937 &generator) {
//   int tries = 1000;

//   qbit listable(qs, N);
//   auto qsRefs = to_ref_wrappers(qs);
//   bool debug = true;

//   for (int k = 0; k < powerOfTwo(N) + tries; k++) {
//     QssIndex idx(N, k);
//     StateVector<N> user_psi = indexToEigen<N>(idx);
//     user_psi = user_psi;
//     if (k >= powerOfTwo(N))
//       user_psi = haar_state_distribution<N>(generator);

//     // Begin prepareState behavior
//     // 1. PrepZ
//     // std::cout << "\nPrepZ is executing\n";
//     eval_hold(qexpr::map(qexpr::_PrepZ, qs));

//     // check simulator against theoretical result
//     // (passes trivially; commented out but left for completeness/posterity)
//     // StateVector<N> computational_zero = indexToEigen<N>(QssIndex(N,0));
//     // assertState<N>(iqs_device, qsRefs, computational_zero, true);

//     /*
//     // Reverse steps of prepareState with invert applied to each
//     //    i.e, the equivalent of
//     //    eval_hold(
//     //      return rotationQMUX<N-1>("Z", qs>>1, qs[0], phis)
//     //      + rotationQMUX<N-1>("Y", qs>>1, qs[0], thetas)
//     //      + disentangleState<N-1>(qs>>1, next_psi, debug)
//     //    );
//     */
//     std::cout << "Performing disentangleState with psi =\n" << user_psi << "\n";

//     // 2.
//     // Classical computation of disentangleState<N-1>
//     StateVector<N - 1> next_psi;
//     std::vector<double> thetas(next_psi.rows());
//     std::vector<double> phis(next_psi.rows());

//     for (unsigned int j = 0; j < next_psi.rows(); j++) {
//       auto [r_c, t_c, theta_c, phi_c] =
//           braketToSpherical(user_psi(2 * j), user_psi((2 * j) + 1), debug);

//       C psi_entry = r_c * std::exp(-i * t_c / 2.0);
//       next_psi(j) = psi_entry;
//       thetas[j] = theta_c;
//       phis[j] = phi_c;
//     }
//     // End:  Classical computation of disentangleState<N-1>

//     // EVALUATE AND CHECK #1:
//     // last term of return statement of disentangleState
//     std::cout << "\ndisentangleState<" << N - 1 << "> is executing\n";
//     // std::cout << "Next psi\n" << next_psi << std::endl;
//     eval_hold(invert(disentangleState<N - 1>(qs >> 1, next_psi, debug)));

//     // theoretical result to check: |next_psi> otimes |0>
//     StateVector<1> zero = indexToEigen<1>(QssIndex(1, 0));
//     StateVector<N> theoretical = kron_vec<N - 1, 1>(next_psi, zero);
//     assertState<N>(iqs_device, qsRefs, theoretical, true);

//     // EVALUATE AND CHECK #2:
//     // middle term of return statement of disentangleState
//     std::cout << "\nrotationQMUX<" << N - 1 << "> Y is executing\n";
//     eval_hold(invert(rotationQMUX<N - 1>("Y", qs >> 1, qs[0], thetas, true)));

//     // theoretical:  YrotationQMUX * |theoretical>
//     theoretical = RYMUXMatrix<N - 1>(thetas).adjoint() * theoretical;
//     assertState<N>(iqs_device, qsRefs, theoretical, true);

//     // EVALUATE AND CHECK #3:
//     // first term of return statement of dientangleState
//     std::cout << "\nrotationQMUX<" << N - 1 << "> Z is executing\n";
//     eval_hold(invert(rotationQMUX<N - 1>("Z", qs >> 1, qs[0], phis, true)));

//     // theoretical:  ZrotationQMUX * |theoretical>
//     theoretical = RZMUXMatrix<N - 1>(phis).adjoint() * theoretical;
//     assertState<N>(iqs_device, qsRefs, theoretical, true);

//     std::cout << "compare running theoretical result to user's inputted psi"
//               << std::endl;
//     bool theo_is_users = isEqualUpToPhase<N>(user_psi, theoretical);
//     if (theo_is_users == false) {
//       std::cout << "The theoretical expression doesn't match the user input.\n";
//       std::cout << "User requested:\n" << user_psi << std::endl;
//       std::cout << "prepareState expressions create:\n"
//                 << theoretical << std::endl;
//       assert(theo_is_users);
//     }

//     std::cout << "compare IQS result to user's inputted psi" << std::endl;
//     assertState<N>(iqs_device, qsRefs, user_psi, true);
//   } // end loop over many states

//   return true;
// }


/// @cond
int main() {
  IqsConfig iqs_config(7, "noiseless", false);
  FullStateSimulator iqs_device(iqs_config);
  QRT_ERROR_T status = iqs_device.ready();
  assert(status == QRT_ERROR_SUCCESS);

  // std::random_device rd;     // random seed
  // std::mt19937 engine(rd()); // random seed for non-deterministic sequence
  std::mt19937 engine(1235); // fixed seed for deterministic sequence

  bool allTestsGood = true;

  /* State Preparation Testing */

  std::cout << "Testing rotationQMUX\n"; // also tested in testPrepareState
  allTestsGood = allTestsGood &&
                 testRotationQMUX<1>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<2>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<3>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<4>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<5>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<6>(iqs_device, engine, 1e2, false);

  std::cout << "Testing disentangleState\n"; // also tested in testPrepareState
  allTestsGood =
      allTestsGood && testDisentangleState<1>(iqs_device, false);
  allTestsGood =
      allTestsGood && testDisentangleState<2>(iqs_device, false);
  allTestsGood =
      allTestsGood && testDisentangleState<3>(iqs_device, false);
  allTestsGood =
      allTestsGood && testDisentangleState<4>(iqs_device, false);
  allTestsGood =
      allTestsGood && testDisentangleState<5>(iqs_device, false);
  allTestsGood =
      allTestsGood && testDisentangleState<6>(iqs_device, false);
  allTestsGood =
      allTestsGood && testDisentangleState<7>(iqs_device, false);

  std::cout << "Testing prepareState\n";
  allTestsGood = allTestsGood &&
                 testPrepareState<1>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testPrepareState<2>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testPrepareState<3>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testPrepareState<4>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testPrepareState<5>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testPrepareState<6>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testPrepareState<7>(iqs_device, engine, 1e2, false);

  if (allTestsGood) {
    std::cout << "All tests passed." << std::endl;
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
/// @endcond
