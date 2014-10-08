INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/mflib
  ${CMAKE_BINARY_DIR}/mflib
)

ADD_LIBRARY(mflib
  MFGraph.cpp
  MFGraph_solve.cpp
  MFSolver.cpp
  MFRegionSolver.cpp
  MFCpgSolver.cpp
  MethylRead.cpp
  MFRegionPrinter.cpp
  MFCpgEstimator.cpp
)

