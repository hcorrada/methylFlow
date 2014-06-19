INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/methylFlow
  ${CMAKE_BINARY_DIR}/methylFlow
)

## Here we define an executable target. Its name is 'lemon-project' and
## is compiled from 'main.cc'. You can add more source files separated
## with whitespaces (including newlines). If you want to build more
## executables, simple repeat (and edit) the following ADD_EXECUTABLE and
## TARGET_LINK_LIBRARIES statements.

ADD_EXECUTABLE(methylFlow
        main.cpp
)

## This tells cmake to install 'lemon-project' to $PREFIX/bin when
## 'make install' is executed. You can give more targets separated
## by whitespaces.

TARGET_LINK_LIBRARIES(methylFlow
  mflib ${LEMON_LIBRARIES} glpk
)


INSTALL(
  TARGETS methylFlow
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)
