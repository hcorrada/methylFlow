include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/avgEvaluation
  ${CMAKE_BINARY_DIR}/avgEvaluation
)

add_executable(avgEvaluation
  main.cpp
)

TARGET_LINK_LIBRARIES(avgEvaluation 
 mflib ${LEMON_LIBRARIES} glpk
) 

INSTALL(
  TARGETS avgEvaluation
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)


