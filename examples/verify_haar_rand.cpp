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
// qabbl/include/
#include <haar_rand.h>

// C++ standard library
#include <iostream>

// Eigen linear algebra library
#include <Eigen/Core>

bool isUnitary(Eigen::MatrixXcd & matrix) {
  double eps = 1e-14;
  Eigen::MatrixXcd check = matrix * matrix.adjoint();

  return check.isIdentity(eps);
}

int main() {
  //an engine that produces a sequence of pseudo-random values.
  std::mt19937 engine(1234);  // fixed seed for deterministic sequence
  
  //std::random_device rd; // random seed 
  //std::mt19937 engine(rd()); // non-deterministic sequence

  int initial_dim_size = 2;
  for(int n = 1; n < 4; n++) {
    int matrix_size = std::pow(initial_dim_size, n);
    std::cout << "Trying " << matrix_size << " x " << matrix_size << std::endl;
    for(int i = 0; i < 1000000; i++) {
      Eigen::MatrixXcd roll = generateHaarDistributedRandomUnitary(matrix_size,
                                                                   engine);
      if(!isUnitary(roll)) {
        std::cout << "Sample " << i << " : failed unitary check!" << std::endl;
        std::cout << roll << std::endl;
        std::cout << roll * roll.adjoint() << std::endl;
        return 1;
      }

      auto test = convertEigenToSTL(roll);
      double eps = 1e-20;
      for(int row=0; row < matrix_size; row++) {
        for(int column = 0; column < matrix_size; column++) {
          //std::cout << test.at(row).at(column) << std::endl;
          if (test.at(row).at(column).real() - roll(row, column).real() > eps) {
            std::cout << "Sample " << i << " : failed conversion check!"
                      << std::endl;
            std::cout << roll << std::endl;
            return 1;
          }
          if (test.at(row).at(column).imag() - roll(row, column).imag() > eps) {
            std::cout << "Sample " << i << " : failed conversion check!"
                      << std::endl;
            std::cout << roll << std::endl;
            return 1;
          }
        }
      }
    }
  }

  return 0;
}
