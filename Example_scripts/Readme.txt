		     RepExpress Example scripts.

This directory contains example script files from the RepExpress
development environment.  The file basic_params.sh contains the
specifications for all of the processing steps.

1. Mapping, counting and gene location step with main script
map_count_reads.sh uses a parameter file like map_params.sh.

2. Combining multiple samples with main script compare_express.sh uses
a parameter file like combine_params.sh for which the file
genloc_file_list.txt is an example of genloc_name_file.

3. Contrasting expression between different sample sets uses a
parameter file like contrast_params.sh where the files
brain_geneloc_list.txt_combined.txt and
testis_genloc_list_txt_combined.txt are output files from step 2
above.

4. Differential expression analysis with DESeq2 under R-Bioconductor
with the script run_deseq2.sh uses the parameter file deseq2_params.sh
and the two TPM featureCounts output files brain_tpm_list.txt &
testis_tpm_list.txt.



