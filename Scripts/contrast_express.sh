#!/bin/bash
#
# RepExpress script to compare TE expression matrices for two
# different compare_express.sh runs which have previously been run
# on all the two sets of samples
#
# Requires two parameters: the name of the basic run parameters
# used for the mapping and the name of a contrast parameter file
# defining the information needed, particularly, the names of
# the compare_express.sh output matrix files and cutoff values
# for the number of samples and the minimum counts for each
# sample.
#
# Peter Stockwell: Mar-2021
#

# check for parameters:

if [[ -z $1 || -z $2 ]]; then

printf "This script needs two parameters:\n";
printf "  1. Name of basic run parameter file\n";
printf "  2. Name of parameter file for contrast run\n";

exit 1;

fi

# Check that we can read this file

if [[ -f $1 && -f $2 ]]; then

# pick up definitions for the run from parameter files

. "$1";
. "$2";

if [[ -n ${verbose} ]]; then

printf "           RepExpress\n";
printf "  TE and gene RNAseq mapping\n\n";
printf "Contrasting Gene Expression hits for two samples\n\n";

printf "Reading basic parameters from '$1'\n";
printf "Reading contrast run parameters from '$2'\n";

printf "Sample 1: from '${TE_matrix1_file}', sample min=${matrix1_sample_min}%%, hit min=${matrix1_hit_min}\n";
printf "Sample 2: from '${TE_matrix2_file}', sample min=${matrix2_sample_min}%%, hit min=${matrix2_hit_min}\n\n";

fi

# Some sanity checks on file names

if [[ -z "${TE_matrix1_file}" ]]; then

printf "Error: no file name defined for TE_matrix1_file in '$2'\n";

exit 1;

fi

if [[ -z "${TE_matrix2_file}" ]]; then

printf "Error: no file name defined for TE_matrix2_file in '$2'\n";

exit 1;

fi

if [[ ! -r "${TE_matrix1_file}" || ! -r "${TE_matrix2_file}" ]]; then

printf "Can't read either or both of '"${TE_matrix1_file}"' & '"${TE_matrix2_file}"\n";

exit 1;

fi

# identify the number of samples for each matrix:

sample1_cols=`awk 'NR==2{print NF-2;}NR>2{exit 0;}' "${TE_matrix1_file}" | head -1`;

sample2_cols=`awk 'NR==2{print NF-2;}NR>2{exit 0;}' "${TE_matrix2_file}" | head -1`;

# take the TE and strand id column from both matrices to make sorted list
# need to have an appropriate header line for later join processing

printf "#TE_id\tstrand\n" > all_TE_list.txt;

cut -f 1,2 "${TE_matrix1_file}" "${TE_matrix2_file}" | sort -u >> all_TE_list.txt;

# awk script to generate the join command

cat << 'JOIN_BLD_CMD' > join_build.awk
# join_build.awk: takes three line file containing:
# 1. Name of sorted TE list.
# 2. Name of TE matrix file 1
# 3. Name of TE matrix file 3
# and generates a join command to
# combine those together to a
# joint matrix.

BEGIN{lcnt = 0;
destfile = "contrast_file.txt";
}

{filelist[lcnt] = $0;
# may as well get the field counts in here.
# could be passed from outside, but this
# should still work OK.
cmd = sprintf("awk '{print NF;}' %s",$0);
cmd | getline fldcount[lcnt];
close(cmd);
lcnt++;
}

END{if (lcnt != 3)
  {
  printf("Error in input: 3 lines needed, got %d\n",lcnt);
  exit 1;
  }
else
  {
  printf("join -o '0,1.2");
  for (col = 3; col <= (fldcount[1] + 0); col++)
    printf(",2.%d",col);
  printf("' %s %s | join -o '0",filelist[0],filelist[1]);
  for (col = 2; col <= fldcount[1]; col++)
    printf(",1.%d",col);
  for (col = 3; col <= (fldcount[2] + 0); col++)
    printf(",2.%d",col);
  printf("' - %s > %s\n",filelist[2],destfile);
  }
}

JOIN_BLD_CMD


# run this script on the files to generate the join script

printf "all_TE_list.txt\n"${TE_matrix1_file}"\n"${TE_matrix2_file}"\n" | \
  awk -f join_build.awk destfile="$2""_joined.txt" > join_matrices.sh;

# run the join script

/bin/sh join_matrices.sh

# awk script to apply limits

cat << 'LIMIT_SCRIPT' > apply_limits.awk
# take values from command line (defaults
# in BEGIN to check that group coverage
# and read hits reach criteria.
#
# Peter Stockwell: Mar-2021

BEGIN{group1_percent = group2_percent = 75.0;
group1_hitmin = group2_hitmin = 50;
samplecnt1 = samplecnt2 = 0;
}

$1~/#/{print $0;}

$1!~/#/{if ((samplecnt1 == 0 || samplecnt2 == 0))
  {
  printf("Error: sample counts not set properly\n");
  exit 1;
  }
else
  {
  validcnts1 = validcnts2 = 0;
  for (i = 3; i <= (samplecnt1 + 2); i++)
    if (($i + 0.0) >= group1_hitmin)
      validcnts1++;
  for (i = (samplecnt1 + 3); i <= NF; i++)
    if (($i + 0.0) >= group2_hitmin)
      validcnts2++;
  if ((validcnts1 >= samplecnt1*group1_percent/100.0) &&
        (validcnts2 >= samplecnt2*group2_percent/100.0))
    print $0;
  }
}
LIMIT_SCRIPT

if [[ -n ${verbose} ]]; then

printf "Selecting joined lines reaching criteria"

fi

printf ""${TE_matrix1_file}" samples=${sample1_cols},\n"${TE_matrix2_file}" samples=${sample2_cols}\n";

awk -f apply_limits.awk samplecnt1=$sample1_cols samplecnt2=$sample2_cols \
  group1_hitmin=$matrix1_hit_min group2_hitmin=$matrix2_hit_min \
  group1_percent=$matrix1_sample_min group2_percent=$matrix2_sample_min \
  "$2""_joined.txt" > "$2""_passed.txt";


if [[ -n ${delete_temp_files} ]]; then

rm join_build.awk
rm apply_limits.awk
rm Join_matrices.sh
#rm "$2""_joined.txt"

fi

if [[ -n ${verbose} ]]; then

printf "\nTE expression matrix build completed.\n\n";
printf "Results written to '"$2""_passed.txt"'\n";

fi

exit 0;

else

printf "Can't open one or both parameter files '$1' or '$2'\n";

exit 1;

fi
