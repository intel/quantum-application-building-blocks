//===--- umatrix.h ----------------------------------------------*- C++ -*-===//
//
//===----------------------------------------------------------------------===//
//
// Copyright 2025 Intel Corporation.
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
/// \file umatrix.h
/// \brief A library of types and functions for working with unitary matrices
///        through the Eigen library.
///
/// Introduces two main data structures, both templated by the number of qubits:
///
///   1. @c StateVector<numQubits> - A 2^n x 1 column vector of complex numbers
///         that represents a qubit register of numQubits size.
///   2. @c UMatrix<numQubits> - A 2^n x 2^n complex unitary matrix that
///         represents a multi-qubit instruction or gate that acts on numQubits.
///
/// These data structures require the number of qubits be supplied at
/// compile-time.
///
//===----------------------------------------------------------------------===//

#ifndef UMATRIX_H
#define UMATRIX_H

/// \def EIGEN_STACK_ALLOCATION_LIMIT
/// Force Eigen to use heap allocation for all matrices (no stack allocation)
#ifndef EIGEN_STACK_ALLOCATION_LIMIT
#define EIGEN_STACK_ALLOCATION_LIMIT 0  // Force all allocations to heap
#endif

// C++ standard library
#include <chrono>
#include <random>

// Intel(R) Quantum SDK header files
#include <qrt_indexing.hpp>
#include <quantum_full_state_simulator_backend.h>

// Eigen linear algebra library
#include <Eigen/Dense>

/// @brief Convenience definition for double-precision complex number
using C = std::complex<double>;

/// @brief Constant to use as a literal for the sqrt(-1)
const C i(0.0, 1.0);

//////////////////
// Declarations //
//////////////////

/// @brief Return the n-th power of two.
/// @param n The exponent.
constexpr int powerOfTwo(int n) { return (n <= 0) ? 1 : 2 * powerOfTwo(n - 1); }


/** @addtogroup qdata Quantum Data
 *  Abstractions for representing quantum data.
 *  @{
 */

/// @brief A type of 2^N vector where N=2^{numQubits}.
template <int numQubits>
using StateVector = Eigen::Vector<C, powerOfTwo(numQubits)>;


/// @brief A type of NxN complex unitary matrices, where
/// N=2^{numQubits}.
///
/// To index into a matrix U, use U(i,j)
/// To initialize a matrix, use
/// UMatrix<numQubits> U {
///    { e11, ..., e1N },
///    ...
///    { eN1, ..., eNN }
/// }
template <int numQubits>
using UMatrix = Eigen::Matrix<C, powerOfTwo(numQubits), powerOfTwo(numQubits),
                              Eigen::RowMajor>;
/** @} */ // end of addtogroup qdata


// Eigen operations:
//  matrix multiplication: m1 * m2 or m1 *= m2
//  conjugate transpose: m.adjoint() or m.adjointInPlace()

/// @brief Generate a unitary matrix of size \c 2^n x \c 2^n from a probability
/// distribution given by the Haar measure.
///
/// @details The inputs are the number of qubits and an
/// initialized \c mt19937 random number generator. The return is an Eigen
/// Matrix object with complex valued entries. The implementation follows the
/// calculation in https://arxiv.org/pdf/math-ph/0609050
///
/// The random number generator needs to be declared in a scope or initialised
/// such that it is persistent through all calls.
///
/// @tparam numQubits the desired size of the unitary matrix
/// @param generator a std::mt19937 random number generator
template <int numQubits>
UMatrix<numQubits> haar_distribution(std::mt19937 &generator);


/// @brief             Return a QssMap object with values for each state key
///                    given by v's entries.
/// @tparam numQubits  Number of qubits in the register represented by v.
/// @param v           The state vector of values.
/// @param TOL         If an entry of v is less than TOL then it is zero.
template <int numQubits>
iqsdk::QssMap<C> eigenToQssMap(const StateVector<numQubits> &v,
                               double TOL = -1.0);


/// @brief             Return a QssMap object with values given by m.
/// @tparam numQubits  Number of qubits in the register represented by m.
/// @param m           A QssMap whose values are used as coefficients in return.
template <int numQubits>
StateVector<numQubits> qssMapToEigen(const iqsdk::QssMap<C> &m);


/// @brief      Return a STL vector of coefficients using the StateVector input.
/// @tparam dim The size of the StateVector.
/// @param v    StateVector to convert to STL vector.
template <typename T, int dim>
std::vector<T> eigenToVector(const Eigen::Vector<T, dim> &v);


/// @brief       Return a StateVector of coefficients using the vector input.
/// @tparam dim  The size of the STL vector.
/// @param v     STL vector to convert to StateVector.
template <typename T, int dim>
Eigen::Vector<T, dim> vectorToEigen(const std::vector<T> &v);


/// @brief Print the probabilities of phi to std out.
/// @tparam numQubits The number of qubits in the register represented by phi.
/// @param phi The StateVector object to print data from.
template <int numQubits> void displayProbabilities(
    const StateVector<numQubits> &phi) {
  auto amps = eigenToQssMap<numQubits>(phi, 0.0001);
  auto probs = iqsdk::FullStateSimulator::amplitudesToProbabilities(amps);
  iqsdk::FullStateSimulator::displayProbabilities(probs);
}


/// \def DEFAULT_EPS
/// Double precision value to use for comparisons; if the difference is smaller
/// than DEFAULT_EPS then the two values are equal.
#define DEFAULT_EPS 0.00001


/// @brief      Check if two complex numbers are equal up to a certain
///             tolerance.
/// @param x    First value for comparison.
/// @param y    Second value for comparison.
/// @param eps  Tolerance to use for deciding if x and y are equal.
bool isEqual(const C &x, const C &y, const double &eps = DEFAULT_EPS);


/// @brief             Check if two state vectors are equivalent up to phase,
///                    given a certain tolerance
/// @tparam numQubits  Number of qubits in the register represented by phi1 and
///                    phi2
/// @param phi1        First state for comparison.
/// @param phi2        Second state for comparison.
/// @param eps         Tolerance for comparing the magnitude of the phase diff
///                    between two states to 1.
template <int numQubits>
bool isEqualUpToPhase(const StateVector<numQubits> &phi1,
                      const StateVector<numQubits> &phi2,
                      const double &eps = DEFAULT_EPS);


/// @brief             Return if the two states phi1 and phi2 the same element
///                    by element.
/// @tparam numQubits  The number of qubits in the register the phi's represent.
/// @param phi1        First state for comparison.
/// @param phi2        Second state for comparison.
template <int numQubits>
bool isEqual(const StateVector<numQubits> &phi1,
             const StateVector<numQubits> &phi2);


/// @brief             Check if two unitary matrices are equal up to a certain
///                    tolerance
/// @tparam numQubits  The number of qubits the matrices operate on.
/// @param mat1        First matrix for comparison.
/// @param mat2        Second matrix for comparison.
/// @param eps         Tolerance to use for deciding if the matrices are equal.
template <int numQubits>
bool isEqual(const UMatrix<numQubits> &mat1, const UMatrix<numQubits> &mat2,
             const double &eps = DEFAULT_EPS);


/// @brief             Compare the multi-qubit operation U to the Identify
///                    operator, i.e. idle or no transformation.
/// @tparam numQubits  Number of qubits the operator U acts on.
/// @param U           Matrix representation of the operator to compare.
/// @param eps         Tolerance for comparing the operator to the Identity
///                    operator.
template <int numQubits>
bool isIdentity(const UMatrix<numQubits> &U, const double &eps = DEFAULT_EPS);


/// @brief      Return true if and only if the input matrix is unitary,
///             e.g. satisfies * Uinv = Udag
/// @tparam     The number of qubits that V operates on.
/// @param V    The multi-qubit instruction to be checked.
/// @param eps  Tolerance to use comparing two numbers; if the difference is
///             smaller than eps the numbers are equal.
template <int numQubits>
bool isUnitary(const UMatrix<numQubits> &V, const double &eps = DEFAULT_EPS);


/// @brief Produce a state vector corresponding to the basis state referred to
/// by the QssIndex i
///
/// @details For example, indexToEigen(QssIndex("|0>")) will produce the state
/// vector {1,0}
template <int numQubits> StateVector<numQubits> indexToEigen(iqsdk::QssIndex i);


/// @brief        Return a UMatrix with elements of Rx with the given rotation
///               angle.
/// @param theta  The size of the rotation.
UMatrix<1> RXMatrix(double theta);


/// @brief        Return a UMatrix with elements of Ry with the given rotation
///               angle.
/// @param theta  The size of the rotation.
UMatrix<1> RYMatrix(double theta);


/// @brief        Return a UMatrix with elements of Rz with the given rotation
///               angle.
/// @param theta  The size of the rotation.
UMatrix<1> RZMatrix(double theta);


/// @brief  UMatrix object with elements of a SWAP gate.
const UMatrix<2> SWAPMatrix {
      { 1, 0, 0, 0 },
      { 0, 0, 1, 0 },
      { 0, 1, 0, 0 },
      { 0, 0, 0, 1 }
};


/// @brief  UMatrix object with elements of a CNOT gate.
const UMatrix<2> CNOTMatrix {
      { 1, 0, 0, 0 },
      { 0, 1, 0, 0 },
      { 0, 0, 0, 1 },
      { 0, 0, 1, 0 }
};


/// @brief  UMatrix object with elements of a Hadamard gate.
const UMatrix<1> HMatrix {
  { 1.0 / std::sqrt(2),   1.0 / std::sqrt(2) },
  { 1.0 / std::sqrt(2),  -1.0 / std::sqrt(2) }
};



/////////////////////
// Implementations //
/////////////////////

template <int numQubits>
UMatrix<numQubits> haar_distribution(std::mt19937 &generator) {
  // 1. Generate an N × N complex matrix Z whose entries are complex standard
  //    normal random variables, where N=2^numQubits.
  std::normal_distribution<double> distribution(1.0, 1.0);
  int N = powerOfTwo(numQubits);

  UMatrix<numQubits> big_z;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      big_z(i, j) = std::complex<double>(distribution(generator) / sqrt(2),
                                         distribution(generator) / sqrt(2));
    }
  }

  // 2. Compute Q and R matrices from Z via linear algebra QR algorithm.
  //    For any complex and square matrix Z, there is a factorization
  //    Z = Q * R where Q is a unitary matrix and R is an uppper triangular
  //    matrix.

  // Expect: Using Ref<> triggers "in-place" algorithms from Eigen to minimize
  // memory footprint throughout the computation.
  Eigen::HouseholderQR<Eigen::Ref<UMatrix<numQubits>>> qr(big_z);
  UMatrix<numQubits> Q = qr.householderQ();
  UMatrix<numQubits> R = qr.matrixQR().template triangularView<Eigen::Upper>();

  // 3. Create the diagonal matrix Λ from Eqn 5.12 where the entries Λ_jj
  //    are the r_jj diagonal elements of R.
  UMatrix<numQubits> lambda;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) {
        lambda(i, j) = R(i, j) / std::abs(R(i, j));
      } else {
        lambda(i, j) = std::complex<double>(0.0, 0.0);
      }
    }
  }

  // 4. The diagonal elements of R′ = Λ^(−1) * R are always real and strictly
  //    positive, therefore the matrix Q′ = QΛ is distributed with Haar measure.

  return Q * lambda;
}

/// Generate a random state vector
template <int numQubits>
StateVector<numQubits> haar_state_distribution(std::mt19937 &generator) {
  UMatrix<numQubits> U = haar_distribution<numQubits>(generator);
  StateVector<numQubits> computational_zero = indexToEigen<numQubits>(iqsdk::QssIndex(numQubits, 0));
  return U * computational_zero;
}

template <int numQubits>
iqsdk::QssMap<C> eigenToQssMap(const StateVector<numQubits> &v, double TOL) {
  iqsdk::QssMap<C> result;
  for (unsigned int i = 0; i < powerOfTwo(numQubits); i++) {
    if (std::sqrt(std::norm(v(i))) >= TOL) {
      result[iqsdk::QssIndex(numQubits, i)] = v(i);
    }
  }
  return result;
}

template <int numQubits>
StateVector<numQubits> qssMapToEigen(const iqsdk::QssMap<C> &m) {
  StateVector<numQubits> result = StateVector<numQubits>::Zero();
  for (auto [key, value] : m) {
    result(key.getIndex()) = value;
  }
  return result;
}

template <typename T, int dim>
std::vector<T> eigenToVector(const Eigen::Vector<T, dim> &v) {
  std::vector<T> result;
  for (unsigned int i = 0; i < dim; i++) {
    result.push_back(v(i));
  }
  return result;
}

template <typename T, int dim>
Eigen::Vector<T, dim> vectorToEigen(const std::vector<T> &v) {
  Eigen::Vector<T, dim> result;
  for (unsigned int i = 0; i < dim; i++) {
    result(i) = v[i];
  }
  return result;
}

template <int numQubits>
StateVector<numQubits> indexToEigen(iqsdk::QssIndex i) {
  iqsdk::QssMap<C> map;
  map[i] = 1.0;
  return qssMapToEigen<numQubits>(map);
}

UMatrix<1> RXMatrix(double theta) {
  UMatrix<1> RX{
    {std::cos(theta / 2.0)    , -i * std::sin(theta / 2.0)},
    {i * std::sin(theta / 2.0),  std::cos(theta / 2.0)}
  };
  return RX;
}
UMatrix<1> RYMatrix(double theta) {
  UMatrix<1> RY{
    {std::cos(theta / 2.0), -std::sin(theta / 2.0)},
    {std::sin(theta / 2.0),  std::cos(theta / 2.0)}
  };
  return RY;
}
UMatrix<1> RZMatrix(double theta) {
  UMatrix<1> RZ{
    {std::exp(-i * theta / 2.0), 0},
    {0                         , std::exp(i * theta / 2.0)}
  };
  return RZ;
}


/// @brief     Return a UMatrix that is the o_plus of U1 and U2 as defined in
///            Shende, V. V., S. S. Bullock, and I. L. Markov. "Synthesis of
///       quantum-logic circuits." IEEE Transactions on Computer-Aided Design of
///            Integrated Circuits and Systems 25.6 (2006): 1000-1010.
/// @tparam numQubits  Number of qubits that U1 and U2 act on.
/// @param U1  One of the "addends" of the direct sum o_plus.
/// @param U2  One of the "addends" of the direct sum o_plus.
template <int numQubits> UMatrix<numQubits + 1> directSum(UMatrix<numQubits> U1, UMatrix<numQubits> U2) {
  constexpr int n = powerOfTwo(numQubits);
  UMatrix<numQubits + 1> result;
  result.template topLeftCorner<n, n>() = U1;
  result.template bottomRightCorner<n, n>() = U2;
  result.template topRightCorner<n, n>() = UMatrix<numQubits>::Zero();
  result.template bottomLeftCorner<n, n>() = UMatrix<numQubits>::Zero();
  return result;
}


/// @brief     Return the matrix that is the kronecker product of U1 and U2.
/// @tparam N  Number of qubits that U1 acts on.
/// @tparam M  Number of qubits that U2 acts on.
/// @param U1  First term in the kronecker product.
/// @param U2  Second term in the kronecker product.
template <int N, int M> UMatrix<N + M> kron(UMatrix<N> U1, UMatrix<M> U2) {
  constexpr int dim2 = powerOfTwo(M);

  UMatrix<N + M> result;
  for (auto i = 0; i < U1.rows(); i++) {
    for (auto j = 0; j < U1.cols(); j++) {
      result.template block<dim2, dim2>(i * dim2, j * dim2) = U1(i, j) * U2;
    }
  }
  return result;
}


/// @brief     Return a \c StateVector that is the kronecker product of the
///            inputted \c StateVectors
/// @tparam N  Number of qubits that V1 acts on.
/// @tparam M  Number of qubits that V2 acts on.
/// @param V1  Left-hand multiplicand.
/// @param V2  Right-hand multiplicand.
template <int N, int M>
StateVector<N + M> kron_vec(const StateVector<N> V1, const StateVector<M> V2) {

  StateVector<N + M> result;
  for (auto i = 0; i < V1.rows(); i++) {
    for (auto j = 0; j < V2.rows(); j++) {
      result(i * V2.rows() + j) = V1(i) * V2(j);
    }
  }
  return result;
}


/// @brief      Return a UMatrix with U11, U12, U21, and U22 as the blocks.
/// @tparam numQubits   Number of qubits U11, U12, U21, and U22 act on.
/// @param U11  The top left block of the output matrix.
/// @param U12  The top right block of the output matrix.
/// @param U21  The bottom left block of the output matrix.
/// @param U22  The bottom right block of the output matrix.
template <int numQubits>
UMatrix<numQubits + 1> blockMatrix(UMatrix<numQubits> U11, UMatrix<numQubits> U12, UMatrix<numQubits> U21,
                           UMatrix<numQubits> U22) {
  constexpr int n = powerOfTwo(numQubits);
  UMatrix<numQubits + 1> result;
  result.template topLeftCorner<n, n>() = U11;
  result.template topRightCorner<n, n>() = U12;
  result.template bottomLeftCorner<n, n>() = U21;
  result.template bottomRightCorner<n, n>() = U22;
  return result;
}




/// @brief      Return the reverse of the input QssIndex.
/// @param idx  The QssIndex to reverse.
iqsdk::QssIndex flipEndianness(iqsdk::QssIndex idx) {
  std::string str = idx.toStringWithoutKet();
  std::reverse(str.begin(), str.end());
  return iqsdk::QssIndex(str);
}


/// @brief         Return a pair of UMatrix objects that has diagonal elements
///                that are Cos or Sin functions with thetas as the angle.
/// @tparam numQubits      Number of qubits the UMatrices act on.
/// @param thetas  A vector of angles, one angle for each diagonal element of
///                the output matrices.
template <int numQubits>
std::pair<UMatrix<numQubits>, UMatrix<numQubits>> cosSinMatrices(std::vector<double> thetas) {
  assert(thetas.size() == powerOfTwo(numQubits));
  UMatrix<numQubits> C = UMatrix<numQubits>::Zero();
  UMatrix<numQubits> S = UMatrix<numQubits>::Zero();

  for (unsigned int i = 0; i < C.rows(); i++) {
    unsigned int j = flipEndianness(iqsdk::QssIndex(numQubits, i)).getIndex();
    C(i, i) = std::cos(thetas.at(j) / 2.0);
    S(i, i) = std::sin(thetas.at(j) / 2.0);
  }
  return std::pair(C, S);
}


/// @brief         Return a UMatric that has diagonal blocks of Cos matrices and
///                off-diagonal blocks that are Sin matrices. The vector of
///                angles are used in both Cos and Sin matrices.
/// @tparam numQubits      Number of qubits the UMatrices act on.
/// @param thetas  A vector of angles, one angle for each diagonal element of
///                the Cos or Sin block matrices.
template <int numQubits> UMatrix<numQubits + 1> cosSinMatrix(std::vector<double> thetas) {
  auto [C, S] = cosSinMatrices<numQubits>(thetas);
  return blockMatrix<numQubits>(C, -S, S, C);
}


/// @brief     Return a UMatrix that represents applying one of the two
///            UMatrix operations in Us; or in other words when applied to input
///            |bs>|phi>, the unitary should produce |bs>(Us[bs]|phi>)
/// @tparam numQubits  Number of qubits that the returned UMatrix operates on.
/// @param Us  A vector of UMatrix operations to choose between
template <int numQubits>
UMatrix<numQubits + 1> controlledUMatrix(const std::vector<UMatrix<1>> &Us) {
  assert(Us.size() == powerOfTwo(numQubits));
  if constexpr (numQubits <= 0) {
    return Us[0];
  } else {
    std::vector<UMatrix<1>> Us0;
    std::vector<UMatrix<1>> Us1;
    for (int i = 0; i < Us.size(); i += 2) {
      Us0.push_back(Us.at(i));
      Us1.push_back(Us.at(i + 1));
    }

    UMatrix<numQubits> U0 = controlledUMatrix<numQubits - 1>(Us0);
    UMatrix<numQubits> U1 = controlledUMatrix<numQubits - 1>(Us1);
    return directSum<numQubits>(U0, U1);
  }
}


/// @brief          Return a UMatrix that represents applying Rz with one of the
///                 angles in thetas.
/// @tparam numQubits       The number of qubits that Rz operates on.
/// @param  thetas  A vector of angles.
template <int numQubits> UMatrix<numQubits + 1> RZMUXMatrix(const std::vector<double> &thetas) {
  if constexpr (numQubits <= 0) {
    return RZMatrix(thetas[0]);
  } else {
    std::vector<double> thetas0;
    std::vector<double> thetas1;
    for (int i = 0; i < thetas.size(); i += 2) {
      thetas0.push_back(thetas.at(i));
      thetas1.push_back(thetas.at(i + 1));
    }

    UMatrix<numQubits> U0 = RZMUXMatrix<numQubits - 1>(thetas0);
    UMatrix<numQubits> U1 = RZMUXMatrix<numQubits - 1>(thetas1);
    return directSum<numQubits>(U0, U1);
  }
}


/// @brief          Return a UMatrix that represents applying Ry with one of the
///                 angles in thetas.
/// @tparam numQubits       The number of qubits that Ry operates on.
/// @param  thetas  A vector of angles.
template <int numQubits> UMatrix<numQubits + 1> RYMUXMatrix(const std::vector<double> &thetas) {
  if constexpr (numQubits <= 0) {
    return RYMatrix(thetas[0]);
  } else {
    std::vector<double> thetas0;
    std::vector<double> thetas1;
    for (int i = 0; i < thetas.size(); i += 2) {
      thetas0.push_back(thetas.at(i));
      thetas1.push_back(thetas.at(i + 1));
    }

    UMatrix<numQubits> U0 = RYMUXMatrix<numQubits - 1>(thetas0);
    UMatrix<numQubits> U1 = RYMUXMatrix<numQubits - 1>(thetas1);
    return directSum<numQubits>(U0, U1);
  }
}


/// @brief          Return a UMatrix that represents applying Rx with one of the
///                 angles in thetas.
/// @tparam numQubits       The number of qubits that Rx operates on.
/// @param  thetas  A vector of angles.
template <int numQubits> UMatrix<numQubits + 1> RXMUXMatrix(const std::vector<double> &thetas) {
  if constexpr (numQubits <= 0) {
    return RXMatrix(thetas[0]);
  } else {
    std::vector<double> thetas0;
    std::vector<double> thetas1;
    for (int i = 0; i < thetas.size(); i += 2) {
      thetas0.push_back(thetas.at(i));
      thetas1.push_back(thetas.at(i + 1));
    }

    UMatrix<numQubits> U0 = RXMUXMatrix<numQubits - 1>(thetas0);
    UMatrix<numQubits> U1 = RXMUXMatrix<numQubits - 1>(thetas1);
    return directSum<numQubits>(U0, U1);
  }
}


/// @brief          Fill the elements of U with values such that applying U will
///                 produce the result.
/// @tparam numQubits       The number of qubits that U operates on.
/// @param  U       The output reference UMatrix.
/// @param idx      The QssIndex for the basis state.
/// @param result   The StateVector that U should produced.
template <int numQubits>
void assignQssIndex(UMatrix<numQubits> &U, iqsdk::QssIndex idx, StateVector<numQubits> result) {
  assert(idx.getNumQubits() == numQubits);
  for (int j = 0; j < powerOfTwo(numQubits); j++) {
    U(idx.getIndex(), j) = result(j);
  }
}


bool isEqual(const C &x, const C &y, const double &eps) {
  bool eq = (std::abs(x.real() - y.real()) < eps);
  eq = eq && (std::abs(x.imag() - y.imag()) < eps);

  return eq;
}


template <int numQubits>
bool isEqual(const StateVector<numQubits> &phi1, const StateVector<numQubits> &phi2) {
  bool eq = true;
  for (int i = 0; i < phi1.rows(); i++) {
    eq = eq && isEqual(phi1(i), phi2(i));
  }

  return eq;
}


template <int numQubits>
bool isEqual(const UMatrix<numQubits> &mat1, const UMatrix<numQubits> &mat2, const double &eps) {
  if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols()) {
    return false;
  }
  for (int i = 0; i < mat1.rows(); i++) {
    for (int j = 0; j < mat1.cols(); j++) {
      if (!isEqual(mat1(i, j), mat2(i, j), eps)) {
        return false;
      }
    }
  }
  return true;
}


// check the equivalence of two states according to the magnitude of the global
// phase difference between the states
//   G. F. Viamontes, I. L. Markov, and J. P. Hayes, “Checking equivalence of
//   quantum circuits and states,” in International Conference on Computer-Aided
//   Design, ICCAD, G. G. E. Gielen, Ed., 2007, pp. 69–74.
template <int numQubits>
bool isEqualUpToPhase(const StateVector<numQubits> &phi1,
                      const StateVector<numQubits> &phi2,
                      const double &eps) {
  return std::abs(phi1.dot(phi2)) - 1 < eps;
}


template <int numQubits>
bool isIdentity(const UMatrix<numQubits> &U, const double &eps) {
  return isEqual<numQubits>(U, UMatrix<numQubits>::Identity(), eps);
}


template <int numQubits>
bool isUnitary(const UMatrix<numQubits> &V, const double &eps) {
  UMatrix<numQubits> result = V * V.adjoint();
  bool b = isIdentity<numQubits>(result, eps);
  result = V.adjoint() * V;
  return b && isIdentity<numQubits>(result, eps);
}


#endif // UMATRIX_H
