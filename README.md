# Quantum Application Building Blocks Library

## About

QABBL is a collection of C++ functions to be used as the building blocks for
applications built using the Intel(R) Quantum SDK.

## Getting Started

Include functions from QABBL in your C++ code project by `#include`-ing the desired header and passing the path to the `qabbl/include` folder to the Intel(R) Quantum SDK's compiler with the `-I` flag.

Example program:
```c++
// my_program.cpp
#include "qalgorithm.h"  // qabbl/include/qalgorithm.h

// remaining code
```

Example Intel(R) Quantum SDK compiler invocation:
```bash
/path/to/IQSDK/intel-quantum-compiler -I/path/to/quantum-library/qabbl/include my_program.cpp
```

Some parts of QABBL are implemented with matrix types from the Eigen linear algebra library. After
getting a copy of Eigen with `git clone https://gitlab.com/libeigen/eigen.git`
you can provide its directory to the Intel(R) Quantum SDK compiler via the `-I` flag.

## Get Help

https://community.intel.com/t5/Intel-Quantum-SDK/bd-p/intel-quantum-sdk

## Essential Documentation

https://developer.intel.com/quantumsdk
