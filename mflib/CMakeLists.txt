INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/mflib
  ${CMAKE_BINARY_DIR}/mflib
)

ADD_LIBRARY(mflib
  MFGraph.cpp
  MFGraph.hpp
  MFGraph_solve.cpp
  MFSolver.cpp
  MFSolver.hpp
  MFRegionSolver.cpp
  MFRegionSolver.hpp
  MFCpgSolver.cpp
  MFCpgSolver.hpp
  MethylRead.cpp
  MethylRead.hpp
  MFRegionPrinter.cpp
  MFRegionPrinter.hpp
  MFCpgEstimator.cpp
  MFCpgEstimator.hpp
)

INSTALL(
  TARGETS mflib
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  COMPONENT library
)		

INSTALL(
  FILES MethylRead.hpp
  DESTINATION ${INSTALL_INCLUDE_DIR}
)

