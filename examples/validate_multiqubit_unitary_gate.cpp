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
/// \file validate_multiqubit_unitary_gate.cpp
/// \brief Demonstration and test of qsd<N>().
///
/// See "multiqubit_gate.h" for detailed descriptions.
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
#include <multiqubit_gate.h>
#include <umatrix.h>
#include <validation_tools.h>


using namespace qlist;
using namespace iqsdk;
using namespace qexpr;


// Check that
// (1) diagonalMUX applied to a qubit register produces the expected state
template <int N>
bool testDiagonalMUX(FullStateSimulator &iqs_device,
                     std::mt19937 &generator, int num_tries = 1,
                     bool debug = true) {
  qbit listable(qs, N + 1);

  // 1. Generate a random diagonal operator D
  UMatrix<N> U = haar_distribution<N>(generator);
  auto [_, D] = diagonalize<N>(U, debug);

  // 2. Calculate E = D \oplus D_dag
  UMatrix<N + 1> E;
  constexpr int n = powerOfTwo(N);
  E.template topLeftCorner<n, n>() = D;
  E.template bottomRightCorner<n, n>() = D.adjoint();
  E.template topRightCorner<n, n>() = UMatrix<N>::Zero();
  E.template bottomLeftCorner<n, n>() = UMatrix<N>::Zero();

  bool b =
      implementsUMatrixHaar<N + 1>(iqs_device, generator, num_tries, qs,
                                   diagonalMUX<N>(D, qs << 1, qs[N]), E, debug);

  if (!b) {
    std::cout << "Failed testDiagonalMUX\n";
    assert(false);
  }
  std::cout << "Passed testDiagonalMUX\n";
  return true;
}


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


template <int N>
bool test_qsd_directSum(FullStateSimulator &iqs_device,
                        std::mt19937 &generator, int num_tries = 1,
                        bool debug = true) {
  qbit listable(qs, N);
  qbit q;

  UMatrix<N> A = haar_distribution<N>(generator);
  UMatrix<N> B = haar_distribution<N>(generator);
  UMatrix<N + 1> V = directSum<N>(A, B);

  QssIndex idx(N + 1, 0);
  auto qsRefs = to_ref_wrappers(qs + QList(q));

  eval_hold(prepareQssIndex<N + 1>(qs + QList(q), idx));
  StateVector<N + 1> sim = indexToEigen<N + 1>(idx);
  assertState<N + 1>(iqs_device, qsRefs, sim);

  // A ⊕ B = (md.U2 ⊕ md.U2) (md.D ⊕ md.D.adjoint()) (md.U1 ⊕ md.U1)
  MUXDecomp<N> md = decomposeMUX<N>(A, B, debug);
  std::cout << "Got decomposition:\n";
  prettyMUXDecomp<N>(md);

  std::cout << "STEP 1\n";
  eval_hold(qsd<N>(qs, md.U1));
  sim = directSum<N>(md.U1, md.U1) * sim;
  assertState<N + 1>(iqs_device, qsRefs, sim);

  std::cout << "STEP 2\n";
  eval_hold(diagonalMUX<N>(md.D, qs, q));
  sim = directSum<N>(md.D, md.D.adjoint()) * sim;
  assertState<N + 1>(iqs_device, qsRefs, sim);

  std::cout << "STEP 3\n";
  eval_hold(qsd<N>(qs, md.U2));
  sim = directSum<N>(md.U2, md.U2) * sim;
  assertState<N + 1>(iqs_device, qsRefs, sim);

  return true;
}


template <int N>
bool testTensorI(FullStateSimulator &iqs_device, std::mt19937 &generator,
                 int num_tries = 1, bool debug = true) {
  qbit listable(qs, N);
  qbit q;

  // Perform H on last qubit

  QssMap<C> psi;
  psi[QssIndex("|00>")] = -0.1965 + 0.577 * i;
  psi[QssIndex("|01>")] = -0.07261 + 0.7028 * i;
  psi[QssIndex("|10>")] = -0.1155 + 0.3393 * i;

  assert(implementsUMatrixOn<N + 1>(
      iqs_device, qs + QList(q), qssMapToEigen<N + 1>(psi),
      _H(q), // SHOULD BE _H(q)
      kron<1, N>(HMatrix, UMatrix<N>::Identity()), true));

  std::cout << "Passed testTensorI\n";
  return true;
}


// Check that
// a pseudo-random unitary matrix is implemented as a
// (1) QSD QExpr applied to each possible basis state of the qubit register
//     produces the expected state
// (2) QSD QExpr applied to a pseudo-random initial state produces the expected
//     states for num_tries generated random initial states.
template <int N>
bool testQSD(FullStateSimulator &iqs_device, std::mt19937 &generator,
             int num_umatrices = 1, int num_tries = 1, bool debug = true) {
  qbit listable(qs, N);
  for (int iter = 0; iter < num_umatrices; iter++) {
    UMatrix<N> U = haar_distribution<N>(generator);
    if (debug) {
      std::cout << "Unitary:\n" << U << "\n";
    }
  // test U on each computational basis vector in the state space
    bool b =
      implementsUMatrix<N>(iqs_device, qs, qsd<N>(qs, U, debug), U, debug);
    if (!b) {
      std::cout << "Failed testQSD (1)\n";
      assert(false);
    }

    // test with randomized state vector
    b = implementsUMatrixHaar<N>(iqs_device, generator, num_tries, qs,
                                 qsd<N>(qs, U, debug), U, debug);
    if (!b) {
      std::cout << "Failed testQSD (2)\n";
      assert(false);
    }
  }
  std::cout << "Passed testQSD for " << num_umatrices << " unitary "
            << N << "-qubit gates, each applied on "
            << powerOfTwo(N) + num_tries << " state vectors\n";
  return true;
}


/* This function is left as historical record for the time being; but it is
   commented out so it doesn't impact compile times.

   Each QExpr was expanded so that the qubit state could be checked between
   every transformation step.
*/
// template <int N>
// bool debugQSDRecursive(FullStateSimulator &iqs_device, UMatrix<N + 1> &U,
//                        StateVector<N + 1> startingState) {
//   bool debug = true;
//   qbit listable(qs, N);
//   qbit q;
//   auto qsRefs = to_ref_wrappers(qs + QList(q));

//   eval_hold(prepareState<N + 1>(qs + QList(q), startingState));
//   StateVector<N + 1> sim = startingState;
//   assertState<N + 1>(iqs_device, qsRefs, sim, true);

//   std::unique_ptr<CSD<N>> decomp = unitaryCSD<N>(U, debug);
//   std::cout << "decomp:\n";
//   prettyCSD<N>(decomp);

//   std::cout << "STEP1: Implementing A2 oplus B2\n\n";
//   eval_hold(qsd_directSum<N>(q, qs, decomp->A2, decomp->B2, debug));
//   sim = (directSum<N>(decomp->A2, decomp->B2)) * sim;
//   assertState<N + 1>(iqs_device, qsRefs, sim, true);

//   /*
//   MUXDecomp<N> md2 = decomposeMUX<N>(decomp->A2, decomp->B2, debug);
//   std::cout << "Decomposed A2 oplus B2 into (md2.U2 oplus md.U2) (md2.D oplus
//   md2.Ddag) (md2.U1 oplus md2.U1):\n"; prettyMUXDecomp<N>(md2); std::cout <<
//   "\n";

//   std::cout << "Implementing md2.U1\n";
//   eval_hold(qsd<N>(qs, md2.U1));
//   sim = (directSum<N>(md2.U1, md2.U1)) * sim;
//   assertState<N+1>(iqs_device, qsRefs, sim);

//   std::cout << "Implementing diagonalMUX with md2.D\n";
//   eval_hold(diagonalMUX<N>(md2.D, qs, q));
//   sim = directSum<N>(md2.D, md2.D.adjoint()) * sim;
//   assertState<N+1>(iqs_device, qsRefs, sim);

//   std::cout << "Implementing md2.U2\n";
//   std::cout << "Current state:\n";
//   displayCurrentProbs(iqs_device, qsRefs);

//   eval_hold(qsd<N>(qs, md2.U2));
//   sim = (directSum<N>(md2.U2, md2.U2)) * sim;
//   assertState<N+1>(iqs_device, qsRefs, sim);
//   */

//   std::cout << "STEP2: Implementing Cos-Sin matrix\n\n";

//   eval_hold(rotationQMUX<N>("Y", qs, q, decomp->thetas));
//   sim = cosSinMatrix<N>(decomp->thetas) * sim;
//   assertState<N + 1>(iqs_device, qsRefs, sim);

//   std::cout << "STEP3: Implementing A1 oplus B1\n\n";
//   eval_hold(qsd_directSum<N>(q, qs, decomp->A1, decomp->B1, debug));
//   sim = (directSum<N>(decomp->A1, decomp->B1)) * sim;
//   assertState<N + 1>(iqs_device, qsRefs, sim);

//   /*
//   std::cout << "Current state:\n";
//   displayCurrentProbs(iqs_device, qsRefs);

//   MUXDecomp<N> md1 = decomposeMUX<N>(decomp->A1, decomp->B1, debug);
//   std::cout << "Decomposed A1 oplus B1 into (md1.U2 oplus md1.U2) (md1.D oplus
//   md1.Ddag) (md1.U1 oplus md1.U1):\n"; prettyMUXDecomp<N>(md1); std::cout <<
//   "\n";

//   std::cout << "Implementing md1.U1\n";
//   eval_hold(qsd<N>(qs, md1.U1));
//   sim = (directSum<N>(md1.U1, md1.U1)) * sim;
//   assertState<N+1>(iqs_device, qsRefs, sim);

//   std::cout << "Executing diagonalMUX with md1.D\n";
//   eval_hold(diagonalMUX<N>(md1.D, qs, q));
//   sim = directSum<N>(md1.D, md1.D.adjoint()) * sim;
//   assertState<N+1>(iqs_device, qsRefs, sim);

//   std::cout << "Executing md1.U2\n";
//   eval_hold(qsd<N>(qs, md1.U2));
//   sim = (directSum<N>(md1.U2, md1.U2)) * sim;
//   assertState<N+1>(iqs_device, qsRefs, sim);
//   */

//   std::cout << "DONE -- checking target state\n";
//   assertState<N + 1>(iqs_device, qsRefs, U * startingState);
//   return true;
// }


/* This function is left as historical record for the time being; but it is
   commented out so it doesn't impact compile times.

   Each QExpr was expanded so that the qubit state could be checked between
   every transformation step.
*/
// bool debugQSD(FullStateSimulator &iqs_device, std::mt19937 &generator) {
//   const int N = 2;
//   qbit listable(qs, N);

//   UMatrix<N> U = haar_distribution<N>(generator);
//   std::cout << "Unitary:\n" << U << "\n\n";

//   StateVector<N> ghz;
//   ghz(0) = 1.0 / std::sqrt(2.0);
//   ghz(ghz.size() - 1) = 1.0 / std::sqrt(2.0);

//   debugQSDRecursive<1>(iqs_device, U, ghz);
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

  /* Multi-qubit Gate Testing */

  std::cout << "Testing rotationQMUX\n";
  allTestsGood = allTestsGood &&
                 testRotationQMUX<1>(iqs_device, engine, 100, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<2>(iqs_device, engine, 100, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<3>(iqs_device, engine, 100, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<4>(iqs_device, engine, 100, false);
  allTestsGood = allTestsGood &&
                 testRotationQMUX<5>(iqs_device, engine, 100, false);

  std::cout << "Testing diagonalMUX\n";
  allTestsGood = allTestsGood &&
                 testDiagonalMUX<1>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testDiagonalMUX<2>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testDiagonalMUX<3>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testDiagonalMUX<4>(iqs_device, engine, 1e2, false);
  allTestsGood = allTestsGood &&
                 testDiagonalMUX<5>(iqs_device, engine, 1e2, false);

  std::cout << "Testing QSD\n";
  allTestsGood = allTestsGood &&
                 testQSD<1>(iqs_device, engine, 1e2, 1e2, false);
  allTestsGood = allTestsGood &&
                 testQSD<2>(iqs_device, engine, 1e2, 1e2, false);
  allTestsGood = allTestsGood &&
                 testQSD<3>(iqs_device, engine, 1e2, 1e2, false);
  allTestsGood = allTestsGood &&
                 testQSD<4>(iqs_device, engine, 1e2, 1e2, false);
  //allTestsGood = allTestsGood &&
  //               testQSD<5>(iqs_device, engine, 1e2, 1e2, false);

  if (allTestsGood) {
    std::cout << "All tests passed." << std::endl;
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
/// @endcond
