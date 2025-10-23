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
/// \file verify_haar_rand.cpp
/// \brief Demonstration and test of generateHaarDistributedRandomUnitary() .
///
/// The example runs in a loop, creating
/// a sequence of unitary matrices for sizes
///    \f$    2^n    \f$.
/// The matrix is confirmed to be unitary to within
///    \f$    10^{-14}    \f$.
//===----------------------------------------------------------------------===//
// qabbl/include/
#include <umatrix.h>

// C++ standard library
#include <iostream>

// Eigen linear algebra library
#include <Eigen/Core>


/// @cond
int main() {
  // an engine that produces a sequence of pseudo-random values.
  std::mt19937 engine(1234); // fixed seed for deterministic sequence

  // std::random_device rd; // random seed
  // std::mt19937 engine(rd()); // non-deterministic sequence

  const int qubit_size = 4; // increase for larger matrices

  std::cout << "Generating random unitary "
            << powerOfTwo(qubit_size) << " x "
            << powerOfTwo(qubit_size) << " matrices." << std::endl;

  for (int i = 0; i < 1000000; i++) {
    UMatrix<qubit_size> roll = haar_distribution<qubit_size>(engine);
    if (!isUnitary<qubit_size>(roll)) {
      std::cout << "Sample " << i << " : failed unitary check!" << std::endl;
      std::cout << roll << std::endl;
      std::cout << "U * U^dagger = " << std::endl;
      std::cout << roll * roll.adjoint() << std::endl;
      return 1;
    }
  }
  std::cout << "All generated matrices are unitary as expected." << std::endl;

  return 0;
}
/// @endcond
