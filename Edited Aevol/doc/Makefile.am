
# Set default behaviour depending on configure options
if WANT_DOXYGEN
all: doxygen
else
all:
endif


# Man pages
dist_man1_MANS = aevol_create.1 aevol_run.1 aevol_modify.1 aevol_propagate.1
dist_man1_MANS += aevol_misc_lineage.1 aevol_misc_ancstats.1
dist_man1_MANS += aevol_misc_robustness.1 aevol_misc_compute_pop_stats.1
dist_man1_MANS += aevol_misc_view_generation.1 aevol_misc_create_eps.1
dist_man1_MANS += aevol_misc_extract.1 aevol_misc_fixed_mutations.1
dist_man1_MANS += aevol_misc_gene_families.1



# How to build the devel doc
doxygen: ../src/*.cpp ../src/libaevol/*.h ../src/libaevol/*.cpp ../src/post_treatments/*.cpp
	doxygen aevol.conf ; echo "document generated\n" > doxygen


# Cleaning rules
CLEANFILES = index.html

clean-local:
	rm -rf html doxygen


# Extra files to include in releases (dist)
EXTRA_DIST = README.txt aevol.conf
