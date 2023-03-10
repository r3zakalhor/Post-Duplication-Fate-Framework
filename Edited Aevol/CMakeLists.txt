# ============================================================================
# Require minimal version of cmake
# ============================================================================
cmake_minimum_required(VERSION 3.8.2)


# ============================================================================
# Set project name and languages
# ============================================================================
project(aevol VERSION 5.0 LANGUAGES C CXX)
add_definitions(-DVERSION="${PROJECT_VERSION}")

# ============================================================================
# Use C++17
# ============================================================================
set(CMAKE_CXX_STANDARD 17)

set(AUTHORIZED_MTPERIOD 607 1279 2281 4253 11213 19937 44497 86243 132049 216091)

set(with-post-treatments ON CACHE BOOL "Whether to build post-treatments")

set(with-x OFF CACHE BOOL "Whether to enable graphical outputs")
set(with-omp ON CACHE BOOL "Whether to enable OpenMP parallelization")

set(enable-profiling OFF CACHE BOOL "Whether to enable profiling")
set(with-tracing OFF CACHE BOOL "Whether to use tracing")
set(enable-perflog OFF CACHE BOOL "Whether to active performance log")
set(enable-normalized-fitness OFF CACHE BOOL "With this option, the NORMALIZED_FITNESS flag is defined, allowing a different fitness calculation")
set(enable-mtperiod 607 CACHE STRING "Period of the Mersen Twister. Autorized values are : 607 1279 2281 4253 11213 19937 44497 86243 132049 216091")
set(enable-trivialjumps OFF CACHE STRING "When this option is set, a trivial jump algorithm will be used instead of the polynomial-based method")
set(enable-mpi OFF CACHE BOOL "When this option is set, the MPI interface is set")
set(enable-omp-sort "" CACHE STRING "Which sorting algorithm to use for sorting mutant individual before processing them (by default, no sorting)")
set(disable-multi-promoters "" CACHE STRING "Forbid that a gene is translated by multiple RNAs: MAX keep only the RNA with maximum production rate (e) or FAREST keep only the RNA which is the farest from the protein")
set(dna-factory-alg "L2G" CACHE STRING "Which memory allocation algorithm to use for managing the DnaFactory pool of DNAs")

set(with-triangle OFF CACHE BOOL "Whether to enable triangle phenotypic target (else Gaussian)")
set(with-raevol OFF CACHE BOOL "Whether to enable R-Aevol (Gene regulation) version")
set(with-detectclone ON CACHE BOOL "Whether to enable clones and not recompute them")
set(with-floatconcentration OFF CACHE BOOL "Whether to enable the encoding of concentration has float (and not double)")
set(with-perf-traces OFF CACHE BOOL "Whether to activate performance traces of (R-)Aevol")
set(with-indiv-perf-traces OFF CACHE BOOL "Whether to activate performance traces (per individual) of (R-)Aevol")

set(search-type "OLD" CACHE STRING "Which SIMD pattern search version of (R-)Aevol")
set(with-opt-diff-search OFF CACHE BOOL "Activate the optimize version of differential search for motifs (RNA and Genes)")
set(without-stats ON CACHE BOOL "Activate (or not) stats (best individuals) for each generation")

set(with-progeny-stats OFF CACHE BOOL "Activate (or not) stats for progeny (should generate a lot of data)")

set(nb-bp "2" CACHE STRING "How many BPs to be used (2 = binary DNA, 4 = classical DNA)")

if ( ${with-indiv-perf-traces} )
        add_definitions(-DWITH_PERF_TRACES_PER_INDIV)
endif ()


if ( ${with-progeny-stats} )
	add_definitions(-DPROGENY_STATS)
endif ()

if ( ${enable-mpi} )
	add_definitions(-DHAVE_MPI)
endif ()

if ( ${with-x} )
    FIND_PACKAGE(X11 REQUIRED)
    if ( X11_FOUND )
        add_definitions(-D__X11)
        INCLUDE_DIRECTORIES(${X11_INCLUDE_DIR})
        LINK_LIBRARIES(${X11_LIBRARIES})
        message("graphical output enable")
    endif ( X11_FOUND )
else ()
    add_definitions(-D__NO_X)
    message("graphical output disable")
endif ()

if(APPLE)
    # X11 headers
    include_directories(AFTER "/opt/X11/include")

    # OpenMP stuff for Clang
    if(CMAKE_C_COMPILER_ID STREQUAL "Clang" OR
       CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
        # C configuration
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp \
            -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        # C++ configuration
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp \
        -Wno-unused-command-line-argument -I /usr/local/opt/libomp/include")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()

    # Hide advanced macOS configuration variables
    mark_as_advanced(
        CMAKE_EXECUTABLE_FORMAT
        CMAKE_INSTALL_PREFIX
        CMAKE_OSX_ARCHITECTURES
        CMAKE_OSX_DEPLOYMENT_TARGET
        CMAKE_OSX_SYSROOT
        Boost_FILESYSTEM_LIBRARY_RELEASE
        Boost_INCLUDE_DIR
    )
endif()

if ( ${with-omp} )
    find_package(OpenMP REQUIRED)

    if ( OPENMP_FOUND )
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    	if(APPLE)
    	    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} /usr/local/lib/libomp.dylib")
	    endif()
    endif ()
endif ()

if ( ${enable-profiling} )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif ()

if ( ${with-tracing} )
	add_definitions(-D__TRACING__)
endif ()

if ( ${without-stats} )
	add_definitions(-DAEVOL_NO_STATS)
endif ()

if ( ${enable-perflog} )
    add_definitions(-D__PERF_LOG__)
endif ()

if ( ${with-perf-traces} )
    add_definitions(-DWITH_PERF_TRACES)
endif ()

if ( ${with-opt-diff-search} )
    add_definitions(-DWITH_OPTIMIZE_DIFF_SEARCH)
endif ()

if ( ${enable-normalized-fitness} )
    add_definitions(-DNORMALIZED_FITNESS)
endif ()

if ( DEFINED enable-mtperiod )
    list(FIND AUTHORIZED_MTPERIOD ${enable-mtperiod} _index)
    if ( ${_index} GREATER -1 )
        set(MTPERIOD ${enable-mtperiod})
    else ()
        message(FATAL_ERROR "period is not a valid Mersenne Twister period")
    endif ()
endif ()
message("Mersen Twister period is set to ${MTPERIOD}")
add_definitions(-DSFMT_MEXP=${MTPERIOD})

if ( ${enable-trivialjumps} )
    if ( ${enable-trivialjumps} MATCHES "[0-9]+" )
        add_definitions(-DTRIVIAL_METHOD_JUMP_SIZE=${enable-trivialjumps})
    else ()
        add_definitions(-DTRIVIAL_METHOD_JUMP_SIZE=1000)
    endif ()
endif ()

if ( ${with-triangle} )
    add_definitions(-DPHENOTYPIC_TARGET_TRIANGLE)
endif ()

if ( ${with-raevol} )
    add_definitions(-D__REGUL)
endif ()

if ( ${with-detectclone} )
    add_definitions(-D__DETECT_CLONE)
endif ()

if ( ${with-floatconcentration} )
    add_definitions(-D__FLOAT_CONCENTRATION)
endif ()

if ( NOT ${enable-omp-sort} STREQUAL "" )
    if ( ${enable-omp-sort} STREQUAL "LDNA" )
        add_definitions(-D__OMP_LIST_SORT=0)
    elseif( ${enable-omp-sort} STREQUAL "SHUFFLE" )
        add_definitions(-D__OMP_LIST_SORT=1)
    endif ()
endif ()

if ( NOT ${disable-multi-promoters} STREQUAL "" )
    if ( ${disable-multi-promoters} STREQUAL "MAX" )
        add_definitions(-D__MULTI_PROMOTERS=0)
    elseif( ${disable-multi-promoters} STREQUAL "FAREST" )
        add_definitions(-D__MULTI_PROMOTERS=1)
    endif ()
endif ()

if ( ${dna-factory-alg} STREQUAL "FIRST" )
    add_definitions(-D__DNA_FACTORY_ALG=0)
elseif ( ${dna-factory-alg} STREQUAL "FIRSTFIT" )
    add_definitions(-D__DNA_FACTORY_ALG=1)
elseif ( ${dna-factory-alg} STREQUAL "L2G" )
    add_definitions(-D__DNA_FACTORY_ALG=2)
elseif ( ${dna-factory-alg} STREQUAL "ALLOCATE" )
    add_definitions(-D__DNA_FACTORY_ALG=3)
endif ()

if ( ${nb-bp} STREQUAL "2" )
    add_definitions(-DBASE_2=1)
elseif ( ${nb-bp} STREQUAL "4" )
    add_definitions(-DBASE_4=1)
endif ()
# ===========================================================================
# Tell CMake to export compile commands for IDE integration
# ===========================================================================
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



# ===========================================================================
# Look for Boost
# ===========================================================================
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.65.0 COMPONENTS filesystem)

# ============================================================================
# Look for Zlib
# ============================================================================
find_package(ZLIB REQUIRED)


# ============================================================================
# Get GNU standard installation directories (GNUInstallDirs module)
# ============================================================================
include(GNUInstallDirs)


# ============================================================================
# Tell CMake where to look for custom modules
# ============================================================================
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# ============================================================================
# Tell cmake where to put binary files.
# By GNU standards "executable programs that users can run" should go in
# bindir a.k.a ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
# and "executable programs to be run by other programs rather than by users"
# in libexecdir a.k.a ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR}
# ============================================================================
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})


# ============================================================================
# Set build type specific compilation flags
# ============================================================================
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -Wall -DDEBUG -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -DDEBUG -g")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math -pipe")
if (NOT CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fexpensive-optimizations")
endif ()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math -pipe")
if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fexpensive-optimizations")
endif ()

# ============================================================================
# Directories containing aevol headers
# ============================================================================
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/mutation")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/phenotype")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/biochemistry")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/rng")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/world")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/population")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/stats")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/cuda")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/io")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/gui")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/utils")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/7")
include_directories("${CMAKE_SOURCE_DIR}/src/libaevol/SFMT-src-1.4")


# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(src)

# ============================================================================
# Include tests
# ============================================================================
add_subdirectory(test EXCLUDE_FROM_ALL)


# ============================================================================
# Adds the 'dist' target (that will use CPack)
# ============================================================================
#add_custom_target(dist COMMAND ${CMAKE_BUILD_TOOL} package_source)


# ============================================================================
# Add the 'uninstall' target (uses a custom script)
# ============================================================================
configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake)
