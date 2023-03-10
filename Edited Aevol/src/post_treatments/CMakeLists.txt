# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================

# ============================================================================
# libaevol-post is a library that contains classes shared by several
# post-treatments
# ============================================================================

include_directories("libaevol_post")
add_subdirectory(libaevol_post)

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_executable(aevol_post_create_eps aevol_post_create_eps.cpp)
target_link_libraries(aevol_post_create_eps PUBLIC aevol)

add_executable(aevol_post_create_eps_7 aevol_post_create_eps_7.cpp)
target_link_libraries(aevol_post_create_eps_7 PUBLIC aevol)

add_executable(aevol_post_create_csv_7 aevol_post_create_csv_7.cpp)
target_link_libraries(aevol_post_create_csv_7 PUBLIC aevol)

add_executable(aevol_post_extract aevol_post_extract.cpp)
target_link_libraries(aevol_post_extract PUBLIC aevol)

add_executable(aevol_post_ancestor_stats aevol_post_ancestor_stats.cpp)
target_link_libraries(aevol_post_ancestor_stats PUBLIC aevol)

add_executable(aevol_post_ancestor_stats_7 aevol_post_ancestor_stats_7.cpp)
target_link_libraries(aevol_post_ancestor_stats_7 PUBLIC aevol)

add_executable(aevol_post_ancestor_mutagenesis aevol_post_ancestor_mutagenesis.cpp)
target_link_libraries(aevol_post_ancestor_mutagenesis PUBLIC aevol libaevol-post)

add_executable(aevol_post_coalescence aevol_post_coalescence.cpp)
target_link_libraries(aevol_post_coalescence PUBLIC aevol)

add_executable(aevol_post_protein_map aevol_post_protein_map.cpp)
target_link_libraries(aevol_post_protein_map PUBLIC aevol)

add_executable(aevol_post_ancestor_robustness aevol_post_ancestor_robustness.cpp)
target_link_libraries(aevol_post_ancestor_robustness PUBLIC aevol libaevol-post)

if ( ${with-raevol} )
add_executable(aevol_post_ancestor_network aevol_post_ancestor_network.cpp)
target_link_libraries(aevol_post_ancestor_network PUBLIC aevol libaevol-post)

add_executable(aevol_post_ancestor_network_knockout_7 aevol_post_ancestor_network_knockout_7.cpp)
target_link_libraries(aevol_post_ancestor_network_knockout_7 PUBLIC aevol)
endif()

add_executable(aevol_post_lineage aevol_post_lineage.cpp)
target_link_libraries(aevol_post_lineage PUBLIC aevol)

add_executable(aevol_post_muller aevol_post_muller.cpp)
target_link_libraries(aevol_post_muller PUBLIC aevol)

add_executable(aevol_post_mutagenesis aevol_post_mutagenesis.cpp)
target_link_libraries(aevol_post_mutagenesis PUBLIC aevol)

add_executable(aevol_post_soft_lineage aevol_post_soft_lineage.cpp)
target_link_libraries(aevol_post_soft_lineage PUBLIC aevol)

add_executable(aevol_post_robustness aevol_post_robustness.cpp)
target_link_libraries(aevol_post_robustness PUBLIC aevol libaevol-post)

add_executable(aevol_post_test_json aevol_post_test_json.cpp)
target_link_libraries(aevol_post_test_json PUBLIC aevol libaevol-post)

add_executable(aevol_post_robustness_json aevol_post_robustness_json.cpp)
target_link_libraries(aevol_post_robustness_json PUBLIC aevol libaevol-post)

add_executable(aevol_post_lineage_json aevol_post_lineage_json.cpp)
target_link_libraries(aevol_post_lineage_json PUBLIC aevol libaevol-post)

add_executable(aevol_post_neutral_mut_acc aevol_post_neutral_mut_acc.cpp)
target_link_libraries(aevol_post_neutral_mut_acc PUBLIC aevol libaevol-post)

add_executable(aevol_post_change_size_neutral_mut aevol_post_change_size_neutral_mut.cpp)
target_link_libraries(aevol_post_change_size_neutral_mut PUBLIC aevol libaevol-post)

add_executable(aevol_post_neutral_shuffle_genome aevol_post_neutral_shuffle_genome.cpp)
target_link_libraries(aevol_post_neutral_shuffle_genome PUBLIC aevol libaevol-post)

add_executable(aevol_post_robustness_bias aevol_post_robustness_bias.cpp)
target_link_libraries(aevol_post_robustness_bias PUBLIC aevol libaevol-post)

add_executable(aevol_post_extract_genome_structure aevol_post_extract_genome_structure.cpp)
target_link_libraries(aevol_post_extract_genome_structure PUBLIC aevol libaevol-post)


add_executable(aevol_post_testJSON testJSON.cpp)
target_link_libraries(aevol_post_testJSON PUBLIC aevol libaevol-post)

# ============================================================================
# Set behaviour on make install
# ============================================================================
