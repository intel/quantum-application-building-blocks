//===--- haar_rand.h ------------------------------*- C++ -*---------------===//
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
/// \file haar_rand.h
/// \brief A library of random sequences of unitary matrices. 
///
/// A function to generate random unitary matrices distributed according to the
/// Haar measure.
///
//===----------------------------------------------------------------------===//

#ifndef HAAR_RAND_H
#define HAAR_RAND_H

// C++ standard library
#include <random>

// Eigen linear algebra library
#include <Eigen/Dense>

/** 
 * @defgroup haarRand Random Unitary Matrices
 * Tools to create sequences of random unitary matrices.
 * @{
 */


/// @brief Generate a unitary matrix of size \c n x \c n from a probability
///        distribution given by the Haar measure.
///
/// @details The inputs are the desired size of the square matrix n and an
/// initialized \c mt19937 random number generator. The return is an Eigen
/// \c Matrix object with complex valued entries. The implementation follows the
/// calculation in https://arxiv.org/pdf/math-ph/0609050
///
/// The random number generator needs to be declared in a scope or initialised
/// such that it is persistent through all calls.
///
/// @param n the desired size of the unitary matrix
/// @param generator a \c std::mt19937 random number generator
Eigen::MatrixXcd generateHaarDistributedRandomUnitary(int n,
                                                      std::mt19937 & generator)
{
  // 1. Generate an N × N complex matrix Z whose entries are complex standard
  //    normal random variables.
  std::normal_distribution<double> distribution(1.0, 1.0);

  Eigen::MatrixXcd big_z;
  big_z.resize(n, n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      big_z(i, j) = std::complex<double>(
          distribution(generator) / sqrt(2),
          distribution(generator) / sqrt(2));
    }
  }

  // 2. Compute Q and R matrices from Z via linear algebra QR algorithm.
  //    For any complex and square matrix Z, there is a factorization
  //    Z = Q * R where Q is a unitary matrix and R is an uppper triangular
  //    matrix.

  // Expect: Using Ref<> triggers "in-place" algorithms from Eigen to minimize
  // memory footprint throughout the computation.
  Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXcd> > qr(big_z);
  Eigen::MatrixXcd Q = qr.householderQ();
  Eigen::MatrixXcd R = qr.matrixQR().triangularView<Eigen::Upper>();

  // 3. Create the diagonal matrix Λ from Eqn 5.12 where the entries Λ_jj
  //    are the r_jj diagonal elements of R.
  Eigen::MatrixXcd lambda;
  lambda.resize(n, n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      if(i == j) { lambda(i, j) = R(i, j) / std::abs(R(i, j)); }
      else { lambda(i, j) = std::complex<double>(0.0, 0.0); }
    }
  }

  // 4. The diagonal elements of R′ = Λ^(−1) * R are always real and strictly
  //    positive, therefore the matrix Q′ = QΛ is distributed with Haar measure.

  return Q * lambda;
}

/** @} */ // end of haarRand


/// @brief Convert a complex Matrix object from Eigen into a vector of vectors
///
/// @details Expects the elements of the input matrix to be std::complex<double>
///
/// @param U_eigen The input matrix to be converted.
std::vector<std::vector<std::complex<double>>>
convertEigenToSTL(Eigen::MatrixXcd & U_eigen) {
  assert(U_eigen.rows() == U_eigen.cols());
  int n = U_eigen.rows();

  // initialize as all 0's to avoid n^2 memory allocations through push_back
  std::vector<std::vector<std::complex<double>>> U_stl(n,
                                          std::vector<std::complex<double>>(n));
  for(int row=0; row < n; row++) {
    for(int column = 0; column < n; column++) {
      U_stl.at(row).at(column) = U_eigen(row, column);
    }
  }

  return U_stl;
}

/** @addtogroup haarRand
 *  
 *  @{
 */

/// @brief Generate a unitary matrix of size \c n x \c n from a probability
/// distribution given by the Haar measure as a vector of vectors
///
/// @details Use the Eigen routines for the generation of the matrix elements
/// and convert to STL type.
/// @return The first index of the returned vector of vectors
/// is treated as i or the row index, and the second as the j or column index.
/// E.G. A.at(0) would return a vector of the full first row.
///      A.at(0).at(3) would return the 3 element of that row, or A_(0,3)

/// @param n the desired size of the unitary matrix
/// @param generator a std::mt19937 random number generator
std::vector<std::vector<std::complex<double>>>
sampleHaarDistributedUnitary(int n, std::mt19937 & generator) {
  Eigen::MatrixXcd U_nn = generateHaarDistributedRandomUnitary(n, generator);

  return convertEigenToSTL(U_nn);
}

/** @} */ // end of addtogroup haarRand

#endif // HAAR_RAND_H
