# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(libaevol)

if (with-post-treatments)
  add_subdirectory(post_treatments)
endif()

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_executable(aevol_create aevol_create.cpp)
target_link_libraries(aevol_create PUBLIC aevol)

add_executable(aevol_modify aevol_modify.cpp)
target_link_libraries(aevol_modify PUBLIC aevol)

add_executable(aevol_run aevol_run.cpp)
target_link_libraries(aevol_run PUBLIC aevol)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(aevol_propagate aevol_propagate.cpp)
    target_link_libraries(aevol_propagate PUBLIC aevol ${Boost_LIBRARIES})
endif()

# ============================================================================
# Set behaviour on make install
# ============================================================================
install(TARGETS aevol_create DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_create)
install(TARGETS aevol_modify DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_modify)
if(Boost_FOUND)
install(TARGETS aevol_propagate DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_propagate)
endif()
install(TARGETS aevol_run DESTINATION ${CMAKE_INSTALL_BINDIR}/aevol_run)