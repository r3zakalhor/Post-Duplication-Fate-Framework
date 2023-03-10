set(AEVOL_REVISION "unknown")

if(EXISTS ${AEVOL_SOURCE_DIR}/.git)
    execute_process(COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${AEVOL_SOURCE_DIR}
        OUTPUT_VARIABLE AEVOL_REVISION_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    execute_process(COMMAND
        git status -s ${AEVOL_SOURCE_DIR}
        WORKING_DIRECTORY ${AEVOL_SOURCE_DIR}
        OUTPUT_VARIABLE AEVOL_SOURCE_MODIFIED
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT AEVOL_SOURCE_MODIFIED STREQUAL "")
        set(AEVOL_REVISION_HASH "${AEVOL_REVISION_HASH}-dirty")
    endif()

    execute_process(COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${AEVOL_SOURCE_DIR}
        OUTPUT_VARIABLE AEVOL_BRANCH_NAME
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    execute_process(COMMAND git log -n 1 --pretty=%cD #--date=short
        WORKING_DIRECTORY ${AEVOL_SOURCE_DIR}
        OUTPUT_VARIABLE AEVOL_REVISION_DATE
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    # If we wanted to use the compile date instead
    # However, this makes the build non-reproducible
    #string(TIMESTAMP AEVOL_COMPILE_DATE "%Y-%m-%d %H:%M:%S")

    set(AEVOL_REVISION "${AEVOL_BRANCH_NAME}, ${AEVOL_REVISION_HASH}, ${AEVOL_REVISION_DATE}")

endif()

if(NOT "${AEVOL_BINARY_DIR}" STREQUAL "")
    configure_file(${AEVOL_BINARY_DIR}/src/libaevol/aevol_version.cpp.in
        ${AEVOL_BINARY_DIR}/src/libaevol/aevol_version.cpp)
endif()
