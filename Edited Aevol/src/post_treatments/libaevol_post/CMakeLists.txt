# ============================================================================
# Libaevol-post is a library that contains classes shared by several
# post-treatments
# ============================================================================

include_directories(".")

add_library(libaevol-post IndivAnalysis.cpp ../../libaevol/io/IOJson.cpp neutral_mutation_exp.cpp
        neutral_mutation_output.cpp  neutral_mutation_output.h Robustness_bias_output.cpp)
target_link_libraries(libaevol-post PRIVATE aevol)
