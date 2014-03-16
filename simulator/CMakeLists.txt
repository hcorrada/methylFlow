include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/mfSimulate
  ${CMAKE_BINARY_DIR}/mfSimulate
)

add_executable(mfSimulate
  main.cpp
  simulator.cpp
)

INSTALL(
  TARGETS mfSimulate
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)


