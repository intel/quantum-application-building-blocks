//===--- multiqubit_gate.h ------------------------*- C++ -*---------------===//
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
/// \file multiqubit_gate.h
/// \brief A quantum expression for a multi-qubit instruction given as a
///        matrix.
///
///        A quantum expression and its helper quantum expressions and functions
///        for decomposing the matrix representation of a multi-qubit gate into
///        a sequence of one- and two-qubit operations.
///
/// Implementation is based on the paper:
///
///   Shende, V. V., S. S. Bullock, and I. L. Markov. "Synthesis of
///   quantum-logic circuits." IEEE Transactions on Computer-Aided Design of
///   Integrated Circuits and Systems 25.6 (2006): 1000-1010.
///
//   Thm 13: Given a 2^(n+1) x 2^(n+1) unitary U, we wish to produce a ZYZ
//           decomposition into four 2^n x 2^n unitaries along with three
//           sets of rotation
//           angles, each of dimension 2^n, such that U is equal to
//
// ----------- Z-MUX(angles1)-------- Y-MUX(angles2)------ Z-MUX(angles2) -----
//                 |                        |                 |
// --/n/--- U1 --- x--------- U2 ---------- x -------- U3 --- x ---------- U4 -
//
//
//===----------------------------------------------------------------------===//

#ifndef MULTI_Q_GATE_H
#define MULTI_Q_GATE_H

// C++ standard library
#include <cassert>
#include <chrono>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

// Intel(R) Quantum SDK header files
#include <clang/Quantum/qexpr.h>
#include <clang/Quantum/quintrinsics.h>
#include <qexpr_utils.h>
#include <qrt_indexing.hpp>
#include <quantum_full_state_simulator_backend.h>

/// \def MKL_Complex16
///  Use STL complex variables with double precision
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

// qabbl/include/
#include <state_prep.h>
#include <umatrix.h>
#include <unitary_ops.h>


//////////////////
// Declarations //
//////////////////

/// @brief Apply the matrix U as a multi-qubit operator by decomposing U into
///        a sequence of one- and two-qubit operations.
///
/// @tparam numQubits  The number of qubits that U operates on.
/// @param qs          The qubits to apply the quantum operation on.
/// @param U           The multi-qubit unitary operation to apply.
/// @param debug       Perform consistency checks and print debug data.
template <int numQubits>
PROTECT QExpr qsd(qlist::QList qs, const UMatrix<numQubits> &U,
                  bool debug = false);


////////////////////
// Implementation //
////////////////////


//////////////////
// Stage 1: Cosine-Sine Decomposition (CSD) //
//////////////////

/// @brief Results from a Cos-Sin Decomposition (CSD).
// Thm 10: CSD
// Given a unitary matrix of dimension 2^{n+1}, return a 4-tuple of unitaries
// (A1, B1, A2, B2), each of dimension 2^n,
// and a vector of 2^n angles `thetas`
// such that
//
// U = ( A1  0  ) V ( A2  0  )
//     ( 0   B1 )   ( 0   B2 )
//
// where
//
// V  = -- Y-MUX(thetas) --
//               |
/// @tparam numQubits  The number of qubits each block matrix of the
///                    decomposition acts on.
template <int numQubits> struct CSD {
  UMatrix<numQubits> A1; //!< Upper block of the left-side matrix.
  UMatrix<numQubits> B1; //!< Lower block of the left-side matrix.
  std::vector<double> thetas; //!< Angles of the Cos-Sin matrix V.
  UMatrix<numQubits> A2; //!< Upper block of the right-side matrix.
  UMatrix<numQubits> B2; //!< Lower block of the right-side matrix.
};


/// @brief             Check that V is the similarity transform of A by
///                    computing D = V^dagger A V and ensuring that D is
///                    diagonal.
/// @tparam numQubits  The number of qubits that A and V operate on.
/// @param A           The matrix or operator to be transformed.
/// @param V           The matrix that applies the similarity transform.
template <int numQubits>
bool checkDiagonalizingSimilarityTransform(const UMatrix<numQubits> &A,
                                           const UMatrix<numQubits> &V) {
  /*
   * Compute D = V^dagger A V to confirm that the matrix of right eigenvectors
   * is a diagonalizing similarity transform. */
  const double eps = 1e-13; // tolerance for equality to 0.0
  auto D = V.adjoint() * A * V;

  for (int i = 0; i < D.rows(); i++) {
    for (int j = 0; j < D.cols(); j++) {
      if (i == j) {
        if (std::abs(D(i, j)) < eps) { // diagonal is almost zero
          std::cout << "D=Vdagger A V failed diagonal check:\n" << D << "\n";
          return false;
        }
      } else {
        if (std::abs(D(i, j)) > eps) { // off-diagonal is non-zero
          std::cout << "D=Vdagger A V failed off-diagonal check:\n"
                    << D << "\n";
          return false;
        }
      }
    }
  }
  return true;
}


//////////
// Shende matrix methods
//////////

/// @brief        Print the matrices of the Cos-Sin Decomposition.
/// @tparam numQubits     The number of qubits each block matrix of the decomposition
///               acts on.
/// @param decomp A pointer to the CSD class.
template <int numQubits> void prettyCSD(const std::unique_ptr<CSD<numQubits>> &decomp) {
  std::cout << "==Printing CSD==\n";
  std::cout << "A1:\n" << decomp->A1 << "\n";
  std::cout << "B1:\n" << decomp->B1 << "\n";
  std::cout << "\nthetas:\n";
  prettyVector(decomp->thetas);
  std::cout << "A2:\n" << decomp->A2 << "\n";
  std::cout << "B2:\n" << decomp->B2 << "\n";
  std::cout << "==End CSD==\n";
}


/// @brief Return a pointer to a Cos-Sin decomposition of the input matrix.
///
// Thm 10: CSD
// Given a unitary matrix of dimension 2^{n+1}, return a 6-tuple of unitaries
// (A1, B1, C, S, A2, B2),
// each of dimension 2^n, such that
//
// U = ( A1  0  ) ( C  -S ) ( A2  0  )
//     ( 0   B1 ) ( S   C ) ( 0   B2 )
//
// AND C and S are real diagonal matrices such that C^2 + S^2 = I
/// @tparam numQubits Number of qubits the matrices of the decomposition act on.
/// @param U      Matrix to be decomposed.
/// @param debug  Perform checks that each matrix of the decomposition is
///               unitary and that the product of the decompositoion give the
///               input U.
template <int numQubits>
std::unique_ptr<CSD<numQubits>> unitaryCSD(const UMatrix<numQubits + 1> &U,
                                   bool debug = false) {
  // 1. Partition X matrix into 4 n/2 by n/2 block matrices
  //    for input to LAPACK routines.
  //    the CS SVD can be used with partitions of general size
  //    determined by the upper left partition's size, p by q
  const int m = powerOfTwo(numQubits + 1); // number of columns or rows in X
  const int p = powerOfTwo(numQubits);     // number of rows of x11 and x12
  const int q = powerOfTwo(numQubits);     // number of columns in x11 and x21

  // array sizes described according to zuncsd documentation to maintain
  // clarity and consistency, see
  // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-0/orcsd-uncsd.html
  // note:  LAPACK_zuncsd() doesn't support unique_ptr type e.g.
  //        std::unique_ptr<C[]> x11(new C[p * p]);
  C *x11 = new C[p * p];
  Eigen::Map<UMatrix<numQubits>> X11(x11);
  X11 = U.template topLeftCorner<p, p>();

  C *x12 = new C[p * (m - q)];
  Eigen::Map<UMatrix<numQubits>> X12(x12);
  X12 = U.template topRightCorner<p, m - q>();

  C *x21 = new C[(m - p) * q];
  Eigen::Map<UMatrix<numQubits>> X21(x21);
  X21 = U.template bottomLeftCorner<p, m - q>();

  C *x22 = new C[(m - p) * (m - q)];
  Eigen::Map<UMatrix<numQubits>> X22(x22);
  X22 = U.template bottomRightCorner<m - p, m - q>();

  // the choices below for leading dimensions ldx## equal to the matrix' number
  // of columns makes them Row-Major format, see
  // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-0/matrix-layout-for-lapack-routines.html
  int ldx11 = q, ldx12 = m - q, ldx21 = q, ldx22 = m - q;

  // 2. peform Cos-Sin Singular Value Decomposition
  double *theta = new double[p];
  C *u1 = new C[p * p];
  C *u2 = new C[(m - p) * (m - p)];
  C *v1h = new C[q * q];
  C *v2h = new C[(m - q) * (m - q)];
  int ldu1 = p, ldu2 = m - p, ldv1h = q, ldv2h = m - q; // Row-Major, see above
  MKL_INT info;
  info = LAPACKE_zuncsd(LAPACK_ROW_MAJOR, 'Y', 'Y', 'Y', 'Y', 'T', 'Z', m, p, q,
                        x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta,
                        u1, ldu1, u2, ldu2, v1h, ldv1h, v2h, ldv2h);

  delete[] x11;
  delete[] x12;
  delete[] x21;
  delete[] x22;

  if (info > 0) {
    printf("The algorithm computing CS-SVD failed to converge.\n");

    exit(1);
  }

  // 3. package the return struct
  std::unique_ptr<CSD<numQubits>> csd = std::make_unique<CSD<numQubits>>();
  csd->A1 = Eigen::Map<UMatrix<numQubits>>(u1);
  csd->B1 = Eigen::Map<UMatrix<numQubits>>(u2);

  csd->A2 = Eigen::Map<UMatrix<numQubits>>(v1h);
  csd->B2 = Eigen::Map<UMatrix<numQubits>>(v2h);

  // The result of the LAPACKE call gives angles thetas in
  // big-endian order e.g. [theta00, theta01, theta10, theta11]
  // but we need thetas in little-endian order e.g.
  // [theta00, theta10, theta01, theta11]
  // csd->thetas.insert(csd->thetas.end(), &theta[0], &theta[p]);
  csd->thetas.resize(p);
  for (unsigned int i = 0; i < p; i++) {
    unsigned int j = flipEndianness(QssIndex(numQubits, i)).getIndex();
    csd->thetas.at(i) = 2.0 * theta[j];
  }

  delete[] u1;
  delete[] u2;
  delete[] theta;
  delete[] v1h;
  delete[] v2h;

  // 4. Check that the results are well-formed

  if (debug) {
    if (!isUnitary<numQubits>(csd->A1)) {
      std::cout << "A1 isn't unitary:\n" << csd->A1 << "\n";
      exit(1);
    }
    if (!isUnitary<numQubits>(csd->A2)) {
      std::cout << "A2 isn't unitary:\n" << csd->A1 << "\n";
      exit(1);
    }
    if (!isUnitary<numQubits>(csd->B1)) {
      std::cout << "B1 isn't unitary:\n" << csd->B1 << "\n";
      exit(1);
    }
    if (!isUnitary<numQubits>(csd->B2)) {
      std::cout << "B2 isn't unitary:\n" << csd->B2 << "\n";
      exit(1);
    }

    // Check that
    // U = ( A1  0  ) V ( A2  0  )
    //     ( 0   B1 )   ( 0   B2 )
    //
    // where
    //
    // V  = -- Y-MUX(thetas) --  =  (  C      S  )
    //               |              (  -S     C  )
    //      ======== x ========
    //
    // and where
    //
    // C = ( cos(theta_0)   0          ... 0                 )
    //     ( -            cos(theta_1) ... 0                 )
    //     ( ... )
    //     ( 0              0          ... cos(theta_{n-1})  )
    // S = ( sin(theta_0)   0          ... 0                 )
    //     ( 0            sin(theta_1) ... 0                 )
    //     ( ... )
    //     ( 0              0          ... sin(theta_{n-1})  )

    UMatrix<numQubits + 1> RHS = directSum<numQubits>(csd->A1, csd->B1) *
                         cosSinMatrix<numQubits>(csd->thetas) *
                         directSum<numQubits>(csd->A2, csd->B2);
    if (!isEqual<numQubits + 1>(U, RHS)) {
      std::cout << "ERROR: unitaryCSD failed to correctly decompose U\n";
      std::cout << "U:\n" << U << "\n";
      prettyCSD(csd);
      std::cout << "RHS:\n" << RHS << "\n";
      assert(0);
    }
  }

  return csd;
}

/////////////////////////////////////////////////////
// Stage 2: Thm 12 -- Demultiplexing a multiplexor //
/////////////////////////////////////////////////////

/// @brief Return both the diagonalized matrix form of the input matrix and
///        the similarity transform matrix.
///
/// Given a unitary matrix V, return a new unitary U and a diagonal unitary D
/// such that V = U D U†
/// @tparam numQubits     The number of qubits that the matrix operates on.
/// @param V      The matrix to diagonlize.
/// @param debug  Print intermediate data during calculation.
template <int numQubits>
std::pair<UMatrix<numQubits>, UMatrix<numQubits>> diagonalize(
        const UMatrix<numQubits> &V,
        bool debug = false) {
  auto rows = V.rows();
  auto cols = V.cols();
  assert(rows > 0);
  assert(cols == rows);
  assert(isUnitary<numQubits>(V));

  /* LAPACKE_zgeev documentation and examples diagonalize input A
   * A = vr^dagger . D . vr
   */

  /* Locals:
   * LAPACKE APIs express matrices as single arrays.
   * For an array a representing a matrix, LAPACK expects an additional int
   * indicating the number of entries in the "leading dimension."
   * This uses Row major representation, see //
   * https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-0/matrix-layout-for-lapack-routines.html
   * so
   *   lda represent the number of entries in a row of matrix a.
   *   ldvl is the number of entries in a row of vl, the left eigenvectors;
   *       unused, but required by LAPACKE.
   *   ldvr is the number of entries in a row of vr, the right eigenvectors.
   */
  MKL_INT n = V.rows(), lda = V.cols(), ldvl = V.cols(), ldvr = V.cols(), info;

  /* Local arrays */
  C *w = new C[n];
  C *vl = new C[ldvl * n];
  C *vr = new C[ldvr * n];

  // We need a copy of a since lapacke_zgeev overwrites a
  C *a = new C[lda * n];
  Eigen::Map<UMatrix<numQubits>> A(a);
  A = V;

  if (debug) { /* Print input matrix */
    std::cout << "A, input to LAPACK_zgeev:\n" << A << "\n";
  }

  /* Solve eigenproblem */
  // Use LAPACKE to solve the eigen-system with zgeev() and return zgeev's
  //    a, the diagonalized matrix data, written on top of the input array
  //   vr, the "right eigen-vectors", the operator that diagonalizes the input
  info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, /* matrix layout */
                       'N',              /* don't write left eigen vectors */
                       'V',              /* write right eigen vectors */
                       n,                /* order of the input matrix */
                       a,                /* array containing the input matrix */
                       lda,              /* leading dim of input matrix */
                       w,                /* array containing eigenvalues */
                       vl,   /* array containing left eigenvectors */
                       ldvl, /* leading dim of left eigenvectors */
                       vr,   /* array containing right eigenvectors */
                       ldvr /* leading dim of right eigenvectors*/);

  if (info > 0) { /* Check for convergence */
    std::cout << "LAPACKE_zgeev failed to compute eigenvalues." << std::endl;
    std::cout << "Input matrix was:\n" << V << "\n";
    std::cout << "Output matrix was:\n"   << A << "\n";
    std::cout << "Eigenvalues were:\n" << w << std::endl;

    exit(1);
  }

  delete[] w;
  delete[] vl;

  // copy the arrays back to UMatrix type
  UMatrix<numQubits> U = Eigen::Map<UMatrix<numQubits>>(vr); // right eigenvectors
  delete[] vr;
  UMatrix<numQubits> D = Eigen::Map<UMatrix<numQubits>>(a);
  delete[] a;

  if (debug) {
    /* Print output matrix */
    std::cout << "Diagonalized operator D:\n" << D << "\n";
    std::cout << "Right eigenvectors:\n" << U << "\n";
  }

  /* Check: is the matrix of right eigenvectors unitary as one expects */
  assert(isUnitary<numQubits>(U));
  return std::make_pair(U, D);
}


/// @brief Abstraction of the decomposition of a direct sum of qubit
///        operators.
///
template <int numQubits> struct MUXDecomp {
  UMatrix<numQubits> U1; //!< Matrix on the right-hand multiplicand.
  UMatrix<numQubits> D;  //!< Diagonal matrix in the center of multiplication.
  UMatrix<numQubits> U2; //!< Matrix of the left-hand multiplicand.
};


/// @brief        Return the decompsition of (V1 ⊕ V2)
///
/// Next we will figure out how to implement (V1 ⊕ V2) for arbitrary unitaries
/// V1 and V2.
/// In particular, we decompose (V1 ⊕ V2) into (U2 ⊕ U2)(D ⊕ D^†)(U1 ⊕ U1)
/// for D a diagonal gate, as described in the proof of Thm 12 (eqn 16).
///
/// It suffices to find U2 and a diagonal matrix D2 that satisfies
///    V1 V2† = U2 D2 U2†
/// Then we can compute D by taking the square root of the entries of D2
/// and compute U1 as U1=D U2† V2
///
/// See https://math.umd.edu/~hking/Normal.pdf
/// @tparam numQubits     The number of qubits that V1 and V2 operate on.
/// @param V1     The left matrix of the direct sum to be decomposed.
/// @param V2     The right matrix of the direct sum to be decomposed.
/// @param debug  Perform a check that the direct sum of the input matrices is
///               equal to the product of the direct sums, i.e.
///               Check that (V1 ⊕ V2) == (U2 ⊕ U2)(D ⊕ D^†)(U1 ⊕ U1)
template <int numQubits>
MUXDecomp<numQubits> decomposeMUX(const UMatrix<numQubits> &V1, const UMatrix<numQubits> &V2,
                          bool debug = false) {

  MUXDecomp<numQubits> result;

  auto [U2, D2] = diagonalize<numQubits>(V1 * V2.adjoint());
  result.U2 = U2;

  // compute D = sqrt(D2)
  for (int i = 0; i < D2.rows(); i++) {
    for (int j = 0; j < D2.cols(); j++) {
      result.D(i, j) = std::sqrt(D2(i, j));
    }
  }

  // Compute U1 = D U2† V2
  result.U1 = result.D * U2.adjoint() * V2;

  if (debug) {
    // Check that  (V1 ⊕ V2)
    // equals      (U2 ⊕ U2)(D ⊕ D^†)(U1 ⊕ U1)
    UMatrix<numQubits + 1> LHS = directSum<numQubits>(V1, V2);
    UMatrix<numQubits + 1> RHS = directSum<numQubits>(result.U2, result.U2) *
                         directSum<numQubits>(result.D, result.D.adjoint()) *
                         directSum<numQubits>(result.U1, result.U1);
    if (!isEqual<numQubits + 1>(LHS, RHS)) {
      std::cout << "ERROR in decomposeMUX; expect the following to be equal:\n"
                << "V1 oplus V2:\n"
                << LHS << "(U1 oplus U2)(D1 oplus D1_dag)(U1 oplus U1):\n"
                << RHS;
      assert(false);
    }
  }

  return result;
}


//////////////////////////////////////
// Stage 3: Diagonal gates as Z-MUX //
//////////////////////////////////////

/// @brief
///
/// Given a diagonal operator D on n qubits, return the QExpr implementing (D ⊕
/// D†)
///
/// Let D have diagonal entries (d{0},...,d{2^n-1}). I claim that
///
/// D ⊕ D†  =   ----- Z-MUX(thetas) ------------
///                         |
///              ---------- x --------------------
///
/// where RZ(theta_i) = ( d{i}     0     )
///                     (  0       d{i}* )
/// i.e. theta_i = -2 i log(d{i})
/// i.e. d{i} = cos(theta_i / 2) - i*sin(theta_i / 2)
///      real(d{i}) = cos(theta_i / 2)
///      imag(d{i}) = -sin(theta_i / 2)
///      imag / real = -tan(theta_i / 2)
///      2* arctan(imag / real) = theta_i
/// @tparam numQubits  Number of qubits that the matrix D operates on.
/// @param D   Matrix to use in direct sum.
/// @param qs  The qubits that the Z-MUX will act on.
/// @param q   The qubit that acts as the control of the MUX.
template <int numQubits>
PROTECT QExpr diagonalMUX(const UMatrix<numQubits> &D, qlist::QList qs, qbit &q) {
  assert(qs.size() == numQubits);
  assert(numQubits > 0);

  std::vector<double> thetas(powerOfTwo(numQubits), 0.0);
  for (size_t i = 0; i < D.rows(); i++) {
    // double angle = -2.0 * i * std::log(D(i,i));
    double angle = 2.0 * std::atan2(-1.0 * D(i, i).imag(), D(i, i).real());
    thetas[i] = angle;
  }

  return rotationQMUX<numQubits>("Z", qs, q, thetas);
}

//////////////////////////////////////////////
// Stage 4: Thm 13: Putting it all together //
//////////////////////////////////////////////
//
// Recall Thm 13: Given a 2^(n+1) x 2^(n+1) unitary U, we wish to produce a
// decomposition into four 2^n x 2^n unitaries along with three sets of
// rotation angles, each of dimension 2^n, such that U is equal to:
//
// ----------- Z-MUX(angles1)-------- Y-MUX(angles2)------ Z-MUX(angles3) -----
//                 |                        |                 |
// --/n/--- U1 --- x--------- U2 ---------- x -------- U3 --- x ---------- U4 -
//


/// @brief     Print the unique matrices of the MUX decomposition.
/// @tparam numQubits  The number of qubits that the matrices act on.
/// @param md  The MUX decomposition to print.
template <int numQubits> void prettyMUXDecomp(const MUXDecomp<numQubits> &md) {
  std::cout << "==Printing MUXDecomp==\n";
  std::cout << "U1: \n" << md.U1 << "\n";
  std::cout << "D: \n" << md.D << "\n";
  std::cout << "U2: \n" << md.U2 << "\n";
  std::cout << "==End MUXDecomp\n";
}


/// @brief Return the quantum expressions for the direct sum A ⊕ B
///
/// @tparam numQubits     Number of qubits A and B each operate on.
/// @param q      The qubit that acts as a control for the diagonal MUX.
/// @param qs     The qubits to appy A and B on.
/// @param A      The left addend of the direct sum.
/// @param B      The right addend of the direct sum.
/// @param debug  Perform consistency checks and print extra steps to stdout.
template <int numQubits>
PROTECT QExpr qsd_directSum(qbit &q, qlist::QList qs, const UMatrix<numQubits> &A,
                            const UMatrix<numQubits> &B, bool debug) {
  assert(qs.size() == numQubits);
  MUXDecomp<numQubits> md = decomposeMUX<numQubits>(A, B, debug);
  return qsd<numQubits>(qs, md.U1)
    + diagonalMUX<numQubits>(md.D, qs, q)
    + qsd<numQubits>(qs, md.U2);
}

/// @brief        Base recursive case for the quantum Shannon decomposition,
/// @param qs     The qubit to apply the quantum instruction to.
/// @param U      The unitary operation to apply.
/// @param debug  Perform consistency checks and print extra steps to stdout.
template <>
PROTECT QExpr qsd<1>(qlist::QList qs, const UMatrix<1> &U, bool debug) {
  // base case
  return fromUnitary1(U, qs[0], debug);
}


template <int numQubits>
PROTECT QExpr qsd(qlist::QList qs, const UMatrix<numQubits> &U, bool debug) {
  // recursive case
  // U = (A1 ⊕ B1) * cosSinMatrix(thetas) * (A2 ⊕ B2)
  std::unique_ptr<CSD<numQubits - 1>> decomp =
    unitaryCSD<numQubits - 1>(U, debug);
  if (debug) {
    std::cout << "Got CSD decomposition:\n";
    prettyCSD<numQubits - 1>(decomp);
  }

  return qsd_directSum<numQubits - 1>(qs[numQubits - 1], qs << 1,
                                      decomp->A2, decomp->B2, debug)
    + rotationQMUX<numQubits - 1>("Y", qs << 1, qs[numQubits - 1],
                                  decomp->thetas)
    + qsd_directSum<numQubits - 1>(qs[numQubits - 1], qs << 1,
                                   decomp->A1, decomp->B1, debug);
}


#endif
