ADD_EXECUTABLE(testMethyl
  testMethyl.cpp
)

TARGET_LINK_LIBRARIES(testMethyl
  mflib
)

configure_file(sim1.tsv sim1.tsv COPYONLY)
configure_file(sim2.tsv sim2.tsv COPYONLY)

add_test(testMethyl testMethyl)
add_test(sim1 ${CMAKE_BINARY_DIR}/methylFlow/methylFlow -i sim1.tsv -o .)
add_test(sim2 ${CMAKE_BINARY_DIR}/methylFlow/methylFlow -i sim2.tsv -o .)

