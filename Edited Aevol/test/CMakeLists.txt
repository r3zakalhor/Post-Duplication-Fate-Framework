# Google Test headers
add_subdirectory(googletest)

# Libraries to link unit tests against
set(test_libs
  aevol
  gtest_main
)

# List of all the tests
set(TESTS
  SampleTest
  IndividualTest
  JumpingMTTest
  TranscriptionTranslationTest
)

# Create a runner for each unit test. Generates a list named TEST_RUNNERS
foreach (TEST IN LISTS TESTS)

  set(TEST_RUNNER run_${TEST})
  set(TEST_RUNNERS ${TEST_RUNNERS} ${TEST_RUNNER})

  add_executable(${TEST} ${TEST}.cpp ../src/post_treatments/aevol_post_robustness_json.cpp)

  # We need the Aevol headers
  target_include_directories(${TEST}
    PRIVATE # these are executables so we don't need to share the includes
    "${CMAKE_SOURCE_DIR}/src/libaevol"
    "${CMAKE_SOURCE_DIR}/src/libaevol/mutation"
    "${CMAKE_SOURCE_DIR}/src/libaevol/phenotype"
    "${CMAKE_SOURCE_DIR}/src/libaevol/biochemistry"
    "${CMAKE_SOURCE_DIR}/src/libaevol/rng"
    "${CMAKE_SOURCE_DIR}/src/libaevol/world"
    "${CMAKE_SOURCE_DIR}/src/libaevol/population"
    "${CMAKE_SOURCE_DIR}/src/libaevol/stats"
    "${CMAKE_SOURCE_DIR}/src/libaevol/cuda"
    "${CMAKE_SOURCE_DIR}/src/libaevol/io"
    "${CMAKE_SOURCE_DIR}/src/libaevol/gui"
    "${CMAKE_SOURCE_DIR}/src/libaevol/utils"
  )

  set_target_properties(${TEST}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )

  target_link_libraries(${TEST} ${test_libs})

  add_custom_target(${TEST_RUNNER}
    COMMENT "Running test ${TEST}..."
    COMMAND ${CMAKE_BINARY_DIR}/test/${TEST}
    DEPENDS ${TEST}
  )

  set_target_properties(${TEST_RUNNER}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )
endforeach(TEST)

# Copy data files for TranscriptionTranslationTest
add_custom_command(TARGET TranscriptionTranslationTest
  PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/test/TranscriptionTranslationTest_files
    ${CMAKE_BINARY_DIR}/test/TranscriptionTranslationTest_files
)


# Add target check (main target for tests)
add_custom_target(check DEPENDS ${TEST_RUNNERS})
