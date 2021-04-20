# contrast_params.sh: defines RepExpress env variables
# to compare TE expression matrices for two
# different compare_express.sh runs which have previously been run
# on all the two sets of samples
#
# Variables that are consistent across a series of runs
# are already defined in basic_params.sh
# 
# Edit this to reflect the parameters required for each contrast run.

# the output directory for the contrast: leave blank for
# output into the current directory.

contrast_output_dir="";

# We need to define the TE expression files containing
# the count matrices from previous compare_express.sh runs.

TE_matrix1_file="brain_genloc_list.txt_combined.txt";

TE_matrix2_file="testis_genloc_list.txt_combined.txt";

# Now define the minimum proportion (%) of samples needed to
# qualify for each set

matrix1_sample_min="75";

matrix2_sample_min="75";

# And then the minimum number of hit counts required for
# each sample - an integer value

matrix1_hit_min="50";

matrix2_hit_min="50";

# The above could have been defined with 1 value each, but
# separating them increases the flexibility.
