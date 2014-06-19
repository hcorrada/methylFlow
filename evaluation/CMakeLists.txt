include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/mfEvaluation
  ${CMAKE_BINARY_DIR}/mfEvaluation
)

add_executable(mfEvaluation
  main.cpp
  evaluation.cpp
)

TARGET_LINK_LIBRARIES(mfEvaluation 
 mflib ${LEMON_LIBRARIES} glpk
) 

INSTALL(
  TARGETS mfEvaluation
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)


