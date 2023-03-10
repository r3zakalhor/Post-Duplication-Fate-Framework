# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(SFMT-src-1.4)

# ============================================================================
# Add targets add their dependencies
# ============================================================================
set(aevol_sources
  ae_enums.h
  ae_logger.h
  AeTime.cpp
  AeTime.h
  aevol.h
  macros.h

  aevol_version.h

  biochemistry/ae_string.cpp
  biochemistry/ae_string.h
  biochemistry/Alignment.cpp
  biochemistry/Alignment.h
  biochemistry/Codon.cpp
  biochemistry/Codon.h
        7/Dna_7.cpp
        7/Dna_7.h
  biochemistry/Dna.cpp
  biochemistry/Dna.h
  biochemistry/GeneticUnit.cpp
  biochemistry/GeneticUnit.h
  biochemistry/Protein.cpp
  biochemistry/Protein.h
  biochemistry/Rna.cpp
  biochemistry/Rna.h
  biochemistry/VisAVis.cpp
  biochemistry/VisAVis.h

  io/Dump.cpp
  io/Dump.h
  io/IOJson.cpp
  io/IOJson.h
  io/json.hpp
  io/LightTree.cpp
  io/LightTree.h
  io/Logging.cpp
  io/Logging.h
  io/OutputManager.cpp
  io/OutputManager.h
  io/ParameterLine.cpp
  io/ParameterLine.h
  io/ParamLoader.cpp
  io/ParamLoader.h
  io/Tree.cpp
  io/Tree.h
  io/DnaReplicationReport.cpp
  io/DnaReplicationReport.h
  io/ReplicationReport.cpp
  io/ReplicationReport.h
#  io/GeneTree.cpp
#  io/GeneTree.h
#  io/GeneTreeNode.cpp
#  io/GeneTreeNode.h

  mutation/Deletion.cpp
  mutation/Deletion.h
        7/DnaMutator.cpp
        7/DnaMutator.h
  mutation/Duplication.cpp
  mutation/Duplication.h
  mutation/HorizontalTransfer.h
  mutation/InsertionHT.cpp
  mutation/InsertionHT.h
  mutation/Inversion.cpp
  mutation/Inversion.h
  mutation/LocalMutation.h
  mutation/Mutation.cpp
  mutation/Mutation.h
        7/MutationEvent.cpp
        7/MutationEvent.h
  mutation/MutationParams.cpp
  mutation/MutationParams.h
  mutation/PointMutation.cpp
  mutation/PointMutation.h
  mutation/Rearrangement.h
  mutation/ReplacementHT.cpp
  mutation/ReplacementHT.h
  mutation/SmallDeletion.cpp
  mutation/SmallDeletion.h
  mutation/SmallInsertion.cpp
  mutation/SmallInsertion.h
  mutation/Translocation.cpp
  mutation/Translocation.h
#  mutation/GeneMutation.cpp
#  mutation/GeneMutation.h

  # The following file does not need compiling as it currently has no methods
  # phenotype/AbstractFuzzy.cpp
  phenotype/AbstractFuzzy.h
  phenotype/Fuzzy.cpp
  phenotype/Fuzzy.h
  phenotype/FuzzyFactory.cpp
  phenotype/FuzzyFactory.h
  phenotype/Gaussian.h
  phenotype/HybridFuzzy.cpp
  phenotype/HybridFuzzy.h
  phenotype/Phenotype.h
  phenotype/PhenotypicSegment.h
  phenotype/PhenotypicTarget.cpp
  phenotype/PhenotypicTarget.h
  phenotype/PhenotypicTargetHandler.cpp
  phenotype/PhenotypicTargetHandler.h
  phenotype/Point.cpp
  phenotype/Point.h

  population/Individual.cpp
  population/Individual.h
  population/IndividualFactory.cpp
  population/IndividualFactory.h
  population/Selection.cpp
  population/Selection.h
        7/Individual_7.cpp
        7/Individual_7.h
        7/Rna_7.cpp
        7/Rna_7.h
        7/Protein_7.cpp
        7/Protein_7.h
        7/Promoter.cpp
        7/Promoter.h
	7/List_Fuzzy.cpp
	7/List_Fuzzy.h
	7/AbstractFuzzy_7.h
	7/FuzzyFactory_7.cpp
	7/FuzzyFactory_7.h
	7/Discrete_Double_Fuzzy.cpp
	7/Discrete_Double_Fuzzy.h

        rng/JumpingMT.cpp
  rng/JumpingMT.h
  rng/JumpPoly.h
  #rng/MersenneTwister.cpp
  #rng/MersenneTwister.h

  stats/AncestorStats.cpp
  stats/AncestorStats.h
  stats/Metrics.h
  stats/NonCodingMetrics.h
  stats/StatRecord.cpp
  stats/StatRecord.h
        7/Stats_7.cpp
        7/Stats_7.h
        7/ExpManager_7.cpp
        7/ExpManager_7.h
  stats/Stats.cpp
  stats/Stats.h
  utils//ObservableEvent.h
  utils/make_unique.h
  utils/Observable.cpp
  utils/Observable.h
  utils/Observer.h
  utils/Utils.cpp
  utils/Utils.h

  world/ExpManager.cpp
  world/ExpManager.h
  world/ExpSetup.cpp
  world/ExpSetup.h
  world/GridCell.cpp
  world/GridCell.h
  world/Habitat.cpp
  world/Habitat.h
  world/HabitatFactory.cpp
  world/HabitatFactory.h
  world/SaveWorld.cpp
  world/SaveWorld.h
  world/World.cpp
  world/World.h
        7/Abstract_Metadata.cpp
        7/Abstract_Metadata.h
      #   7/Map_Metadata.cpp
      #   7/Map_Metadata.h
      #   7/DynTab_Metadata.cpp
      #   7/DynTab_Metadata.h
        7/List_Metadata.cpp
        7/List_Metadata.h
        7/DnaFactory.cpp
        7/DnaFactory.h
        7/Vector_Fuzzy.h
        7/Vector_Fuzzy.cpp
  )

if ( ${with-x} )
    set (aevol_sources ${aevol_sources}
        gui/X11Window.cpp
        gui/X11Window.h
        population/Individual_X11.cpp
        population/Individual_X11.h
        world/ExpManager_X11.cpp
        world/ExpManager_X11.h
    )
endif ()

if ( ${with-raevol} )
    set(aevol_sources ${aevol_sources}
          raevol/Habitat_R.h
          raevol/Habitat_R.cpp
          raevol/Individual_R.h
          raevol/Individual_R.cpp
          raevol/PhenotypicTarget_R.h
          raevol/PhenotypicTarget_R.cpp
          raevol/PhenotypicTargetHandler_R.h
          raevol/PhenotypicTargetHandler_R.cpp
          raevol/Protein_R.h
          raevol/Protein_R.cpp
          raevol/Rna_R.h
          raevol/Rna_R.cpp
	  #"SIMD_Abstract Metadata.cpp" "SIMD_Abstract Metadata.h" SIMD_Map_Metadata.cpp SIMD_Map_Metadata.h
            raevol/SIMD_PhenotypicTargetHandler_R.cpp raevol/SIMD_PhenotypicTargetHandler_R.h)
    if ( ${with-x})
        set(aevol_sources ${aevol_sources}
            raevol/Individual_R_X11.h
            raevol/Individual_R_X11.cpp)
    endif ()

endif ()

# ============================================================================
# Generate version file aevol_version.cpp
# ============================================================================

configure_file(aevol_version.cpp.in ${CMAKE_BINARY_DIR}/src/libaevol/aevol_version.cpp.in @ONLY)

# configure further at build time (always built as run_always.txt is never actually created)
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/src/libaevol/aevol_version.cpp run_always.txt
    COMMAND ${CMAKE_COMMAND}
        -DAEVOL_DIR=${PROJECT_SOURCE_DIR}
        -DAEVOL_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DAEVOL_BINARY_DIR=${CMAKE_BINARY_DIR}
        -P ${CMAKE_SOURCE_DIR}/CMakeScripts/aevol_version.cmake
    COMMENT "Generating aevol_version.cpp"
)

list(APPEND aevol_sources
    ${CMAKE_BINARY_DIR}/src/libaevol/aevol_version.cpp)

add_library(aevol ${aevol_sources})

# ============================================================================
# Use GSL and zlib
# ============================================================================
target_link_libraries(aevol PUBLIC ${ZLIB_LIBRARY})
target_link_libraries(aevol PUBLIC sfmt)


# ============================================================================
# Make STDC MACROS available (for fixed width integers)
# ============================================================================
target_compile_definitions(aevol
  PUBLIC __STDC_FORMAT_MACROS
  PUBLIC __STDC_CONSTANT_MACROS)
