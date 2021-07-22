#!/bin/bash
#
# RepExpress script to run DESeq2 under R on
# featureCounts output files for unique expression.
# This script runs on _U_FC.tpm files
# generated in the map_count_reads.sh operation, which
# will previously have been run on all samples.
#
# Requires two parameters: the name of the basic run parameters
# used for the mapping and the name of a contrast parameter file
# defining the information needed, particularly, the names of
# files of lists of the data sets to be compared by DESeq2
#
# Rscript work based on https://lashlock.github.io/compbio/R_presentation.html
#
# Peter Stockwell: May-2021
#

# check for parameters:

if [[ -z $1 || -z $2 ]]; then

printf "This script needs two parameters:\n";
printf "  1. Name of basic run parameter file\n";
printf "  2. Name of parameter file for DEseq2 run\n";

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
printf "Running DESeq2 on two sample sets\n\n";

printf "Reading basic parameters from '$1'\n";
printf "Reading DESeq run parameters from '$2'\n";

printf "Sample Set 1: from '${deseq2_file_list_1}' with name '${sample_name_1}'\n";
printf "Sample Set 2: from '${deseq2_file_list_2}' with name '${sample_name_2}'\n";

fi

# Some sanity checks on file names

if [[ -z "${deseq2_file_list_1}" ]]; then

printf "Error: no file name defined for deseq2_file_list_1 in '$2'\n";

exit 1;

fi

if [[ -z "${deseq2_file_list_2}" ]]; then

printf "Error: no file name defined for deseq2_file_list_2 in '$2'\n";

exit 1;

fi

if [[ ! -r "${deseq2_file_list_1}" || ! -r "${deseq2_file_list_2}" ]]; then

printf "Can't read either or both of '"${deseq2_file_list_1}"' & '"${deseq2_file_list_2}"'\n";

exit 1;

fi

# take the TE id column from all files to make unique sorted list

if [[ -n ${verbose} ]]; then

printf "Creating list of unique sorted TE IDs\n This may take some time\n";

fi

awk 'BEGIN{printf("tail -n +3");}{printf(" %s",$0);}END{printf(" | cut -f 1 | sort -u | tail -n +2 > all_deseq2_TEs.txt\n");}' \
  "${deseq2_file_list_1}" "${deseq2_file_list_2}" | /bin/sh

# run this script on the files to generate the join script

# awk script to make the matrix and metadata files

cat << 'MATRIX_BUILD' > build_deseq2_files.awk
# build_deseq2_files.awk: script to take two input files of
# featureCounts output files and generate join commands to
# build the two input files required for DESeq2 running
# under R.

BEGIN{file1_descr = "control";
file2_descr = "treatment";
meta_file = "DESeq2_metadata.txt";
printf("ColID\tdescr\n") > meta_file;
nullfld = "0";
fcount = 0;
id_list_file = "all_deseq2_TEs.txt";
ostring = "0,";
countcol = 7;
basecmd = sprintf("join --nocheck-order -a 1 -e \"%s\"",nullfld);
prvcmd = sprintf("%s -o '%s2.%d' %s ",basecmd,ostring,countcol,id_list_file);
outfile = "DESeq2_count_matrix.txt";
cmd = "";
prvfile = "";
hdr = "#TE_id";
}

FNR==NR{printf("%s\t%s\n",basename($0),file1_descr) > meta_file;
add_join($0);
}

FNR < NR{printf("%s\t%s\n",basename($0),file2_descr) > meta_file;
add_join($0);
}

function basename(fullpath)
{
ns = split(fullpath,fps,"/");
if (ns > 0)
  return(fps[ns]);
else
  return("");
}

function add_join(jfilename)
# to add the next join section to existing
# command string
{
hdr = hdr "\t" basename(jfilename);
ostring = "0,";
for (i = 1; i <= NR; i++)
  ostring = ostring sprintf("1.%d,",i+1);
cmd = cmd prvcmd;
prvcmd = sprintf(" '%s' | %s -o '%s2.%d' - ",jfilename,basecmd,ostring,countcol);
prvfile = jfilename;
}

END{printf("printf \"%s\\n\" > %s\n",hdr,outfile);
printf("%s '%s' | tr \" \" \"\\t\" | awk '$1!~/==>/&&$1!=\"Geneid\"' >> '%s'\n",cmd,prvfile,outfile);
}
MATRIX_BUILD

if [[ -n ${verbose} ]]; then

printf "Creating DESeq2 input files for matrix and metadata\n";

fi

#awk -f build_deseq2_files.awk file1_descr="${sample_name_1}" file2_descr="${sample_name_2}" "${deseq2_file_list_1}" "${deseq2_file_list_2}" | /bin/sh

awk -f build_deseq2_files.awk file1_descr="${sample_name_1}" file2_descr="${sample_name_2}" "${deseq2_file_list_1}" "${deseq2_file_list_2}" > do_the_join.sh

/bin/sh ./do_the_join.sh

if [[ ! -n ${DESeq2_output_file} ]]; then

DESeq2_output_file=$(basename "${deseq2_file_list_1}")"_"$(basename "${deseq2_file_list_2}")".txt";

fi

# easier to generate the R script in two steps in order to avoid losing '$padj'

cat << 'DESEQ_RSCRIPT' > run_DESeq2.R
library("DESeq2")
countData <- read.csv("DESeq2_count_matrix.txt", header=TRUE, sep="\t")
metaData <- read.csv("DESeq2_metadata.txt", header=TRUE, sep="\t")
dds <- DESeqDataSetFromMatrix(countData=countData,colData=metaData,design=~descr,tidy=TRUE)
dds <- DESeq(dds)
res <- results(dds)
res_sorted <- res[order(res$padj),]
DESEQ_RSCRIPT

cat << RSCRIPT_ADDENDUM >> run_DESeq2.R
write.table(res_sorted,file="${DESeq2_output_file}",sep="\t",quote=FALSE,col.names=NA,row.names=TRUE,append=TRUE)
q()
RSCRIPT_ADDENDUM

if [[ -n ${verbose} ]]; then

printf "Running Rscript on 'run_DESeq2.R' to generate DESeq2 output\n";
printf " This may take some minutes to complete\n";

fi

# put col1 header into output file, since R doesn't

printf "TE_id" > "${DESeq2_output_file}";

# run R on the script

Rscript run_DESeq2.R;

if [[ -n ${verbose} ]]; then

printf "DESeq2 run completed, results written to '${DESeq2_output_file}'\n\n";

fi

if [[ -n ${delete_temp_files} ]]; then

rm build_deseq2_files.awk
rm run_DESeq2.R
rm do_the_join.sh

fi

exit 0;

else

printf "Can't open either or both of '$1' or '$2'\n";

exit 1;

fi
