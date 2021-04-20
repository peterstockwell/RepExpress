# basic_params.sh: define basic RepExpress env variables to show:
# 
# 1. Control verbosity as script runs
# 2. Locations of required executables (e.g. STAR, Stringtie, featureCounts)
# 3. Locations of raw genomic sequence file(s) and generated index for mapping:
# 4. Locations of gtf annotation files
#
# Edit this file to reflect the parameters required on your target system.
# This information is expected to suit a whole series of runs which would
# be made against the same genome and annotations.  Definitions for each
# individual run should be made in a run parameter script (run_params.sh,
# for instance).

# 1. Control verbosity: set to empty string to reduce feedback

verbose="yes";

# Control retention of working scripts, set this string to empty to retain them

delete_temp_files="yes";

# 2.  Locations of required executables and invariant run parameters:

# path to STAR executable: can be left empty if STAR is already on your exec PATH
#  The command 'which STAR' will indicate this for you.

path_to_star="";

# STAR mapping run parameters:

starfiltermax="150";
staranchormax="150";
starthreads="4";

# path to featureCounts executable: empty if on your path
# 'which featureCounts' will indicate this

path_to_featurecounts="";

# featureCounts parameters

featurecounts_threads="4";
featurecounts_overlap="--minOverlap 25";

# path to stringtie executable: empty if already on your path
# 'which stringtie' will indicate this

path_to_stringtie="";

# stringtie parameters

stringtie_threads="4";

# path to DMAP executables, particularly identgeneloc.
#  leave empty if these are already on your path
# 'which identgeneloc' will indicate this

path_to_dmap="";

# 3. Locations of raw genomic sequence file(s) and generated index for mapping:

# genome fasta files - this is for all sequences in one large fasta file.
#  It is possible to have separate files for each chromosome, but it
#  tends to get a bit messy.

genome_fasta_file="/mnt/hcs/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/hs_gencode_GRCh37/GRCh37.primary_assembly.genome.fa.gz"

# the location where the index files are written and read from

star_genome_dir="/mnt/hcs/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/sra_data/RepExpress_building/dsm_hs_ref_GRCh37/";

# 4. Locations of gtf annotation files

# name of the UCSC repeat element source file downloaded from
# https://genome.ucsc.edu/cgi-bin/hgTables
# with web page settings:
#
# 'assembly': must be set to the required build
#
# 'group': set to Repeats
#
# 'output format': set to 'all fields from selected table' in order to retrieve
#   details that are omitted from the GTF - gene transfer format (limited) format.
#
# 'output file': should be set to something to avoid having the data sent directly
#  for display by the browser (e.g. hg19_ucsc_repeats.txt).
#
# The file can be gzip compressed or not.  The gzip compression at the web interface
#  didn't seem to work.

ucsc_repeat_src="/mnt/hcs/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/sra_data/STAR_chi/FeatureCount/hg19_ucsc_repeats.txt.gz"

# name of gencode gene annotation gtf file, available from
# https://www.gencodegenes.org/human/release_32lift37.html or
# by ftp from ftp.ebi.ac.uk at /pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping
#
# The desired file for GRCh37/hg19 is gencode.v32lift37.annotation.gtf.gz
#
# The GRCh38 release is available from the related directory.
#

gencode_gene_gtf_src="/mnt/hcs/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/hs_gencode_GRCh37/gencode.v32lift37.annotation.gtf"

# Name of dir to save repeat and gencode gtfs: blank will use current default
#  this assumes both are going to be in the same location.

repeat_gene_gtf_dir="/mnt/hcs/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/sra_data/RepExpress_building/dsm_hs_ref_GRCh37/";

# Name of combined gene+repeat gtf file - leave blank for a default name

gene_repeat_gtf=""

# Stuff below here is combining and processing information from above.
# It should not be necessary to change anything below

# check the gencode gene gtf details:
# if it has a gzip extension then modify the name appropriately

if [[ "${gencode_gene_gtf_src}" == *".gz" ]]; then

  gencode_gene_gtf="${repeat_gene_gtf_dir}"$(basename "${gencode_gene_gtf_src}" ".gz");
  

else

  gencode_gene_gtf="${gencode_gene_gtf_src}";

fi

if [[ -z ${gene_repeat_gtf} ]]; then

gene_repeat_gtf="${repeat_gene_gtf_dir}""GenesPlusRepeats.gtf"

fi

# name for ensembl ID vs gene name file

ensid_vs_gname="${repeat_gene_gtf_dir}""ensid_vs_gname.txt";

# name for unique repeat gtf file:

ucsc_repeats_uniq_gtf="${repeat_gene_gtf_dir}""ucsc_repeats_uniq.gtf"

