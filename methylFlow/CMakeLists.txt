INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
)

## Here we define an executable target. Its name is 'lemon-project' and
## is compiled from 'main.cc'. You can add more source files separated
## with whitespaces (including newlines). If you want to build more
## executables, simple repeat (and edit) the following ADD_EXECUTABLE and
## TARGET_LINK_LIBRARIES statements.

ADD_EXECUTABLE(methylFlow
        main.cpp
	MFGraph.cpp
	MFGraph_solve.cpp
	MethylRead.cpp
	MFRegionPrinter.cpp
)

## This tells cmake to install 'lemon-project' to $PREFIX/bin when
## 'make install' is executed. You can give more targets separated
## by whitespaces.

TARGET_LINK_LIBRARIES(methylFlow
  ${LEMON_LIBRARIES} glpk
)

INSTALL(
  TARGETS methylFlow
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)
