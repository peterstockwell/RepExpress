# map_params.sh: to define RepExpress env variables
# for individual mapping runs
#
# Variables that are consistent across a series of runs
# are already defined in basic_params.sh
# 
# Edit this to reflect the parameters required for each mapping run.

# the output directory for STAR mapping results: leave blank for
# output into the current directory.  Don't put the value "./"
# there, since this will cause the stringtie run to fail.

#mapping_output_dir="./star_out/";
mapping_output_dir="";

# Your data files for mapping: Read1 required, Read 2 optional, leave blank if not used
# can be gzip compressed if your file system allows this.

read1_fastq="/mnt/hcs/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/sra_data/SRR2422921/SRR2422921_1_at.fastq"
read2_fastq="/mnt/hcs/dsm-pathology-ecclesRNA/Erin_Macaulay_placental/sra_data/SRR2422921/SRR2422921_2_at.fastq"

# featureCounts needs a strandedness parameter reflecting the way the
# library was generated.  Values are:
#   unstranded       "-s 0"  (default if setting is left blank, necessary for single ended reads)
#   stranded         "-s 1"
#   reverse stranded "-s 2"

featurecounts_strandedness="-s 1"
