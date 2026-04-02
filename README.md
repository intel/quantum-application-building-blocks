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
/path/to/IQSDK/intel-quantum-compiler -I/path/to/quantum-application-building-blocks/include my_program.cpp
```

Some parts of QABBL are implemented with code from third-party C++ libraries.
These libraries can be found in the list of include directives.
After getting a copy of each library listed in the program (typically through git clone),
compile the program by modifying the compilation command by adding an entry per library of the `-I` flag and the path to the library (in the same way as done for the `qabbl/include` as described above).

## Get Help

https://community.intel.com/t5/Intel-Quantum-SDK/bd-p/intel-quantum-sdk

## Essential Documentation

https://intel.github.io/quantum-sdk-docs/
