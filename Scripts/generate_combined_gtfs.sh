#!/bin/bash
#
# RepExpress genome and gtf building script to perform:
#  1: build GTF files ucsc repeats with unique repeat IDs
#  2: combine with gencode gene information
#  3: build STAR genome index
#
# requires a parameter: the name of a file containing information on
# paths of executables, run time parameters and genome and annotation
# files.
#
# Peter Stockwell: Feb-2021
#

# check for parameter:

if [[ -z $1 ]]; then

echo "This script needs a parameter:";
echo "  1. Name of basic parameter file";

exit 1;

fi

# Check that we can read this file

if [[ -f $1 ]]; then

# pick up definitions from the parameter file

if [[ -n ${verbose} ]]; then

printf "           RepExpress\n";
printf "  TE and gene RNAseq mapping\n\n";
printf "  STAR genome and gtf building phase\n\n";
printf "Reading basic information from '$1'\n";

fi

. "$1";

# make awk scripts for further processing

cat << 'SCRIPT1' > make_uniq_gtf_ex_allfields.awk
# make_uniq_gtf_ex_allfields.awk: script to take UCSC 'all fields'
# output from https://genome.ucsc.edu/cgi-bin/hgTables
# and generate a gtf equivalent with uniq repeat IDS and
# with meaningful class_id values.
#
# this uses an awk array based method of selecting relevant lines.
# 
# also puts gene_biotype value to avoid STAR's Unknown biotype
# geneInfo.tab messages
#
# Peter Stockwell: 19-Mar-2020
# 8-Apr-2020 version includes 'chr' in chromosome IDs and
# rejects fix & random chromosome fragments.
#

BEGIN{okclass["SINE"] = okclass["LINE"] = okclass["LTR"] = okclass["DNA"] = 1;}

$1!~/#/&&$6!~/_/{if ($12 in okclass)
  {
  uniqid = sprintf("%s_%s_%s_%s_%s",substr($6,4),$7,$8,$11,$12);
  printf("%s\tUCSCrepeats\texon\t%s\t%s\t%s.0\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\";",
           $6,$7,$8,$2,$10,uniqid,uniqid);
  printf(" gene_biotype \"repeat\"; family_id \"%s\"; class_id \"%s\";\n",$11,$12);
  }
}

SCRIPT1

# Apply this script to UCSC repeat file, checking if gzip compressed.

if [[ -n ${verbose} ]]; then

printf "Creating repeat gtf file '${ucsc_repeats_uniq_gtf}'\n";
printf "with unique IDs from '${ucsc_repeat_src}'\n";

fi

if [[ "${ucsc_repeat_src}" == *".gz" ]]; then

gzip -dc "${ucsc_repeat_src}" | awk -f make_uniq_gtf_ex_allfields.awk > "${ucsc_repeats_uniq_gtf}";

else

awk -f make_uniq_gtf_ex_allfields.awk "${ucsc_repeat_src}" > "${ucsc_repeats_uniq_gtf}";

fi

# check for uncompressed gencode gene gtf file - if not generate it.

if [[ ! -f "${gencode_gene_gtf}" ]]; then

  gzip -dc "${gencode_gene_gtf_src}" > "${repeat_gene_gtf_dir}""${gencode_gene_gtf}";

fi

# append modified repeat gtf to gencode genomic gtf

cat "${gencode_gene_gtf}" "${ucsc_repeats_uniq_gtf}" > "${gene_repeat_gtf}";

if [[ -n ${verbose} ]]; then

printf "Files '${gencode_gene_gtf}' and\n";
printf " '${ucsc_repeat_src}' combined to\n";
printf " '${gene_repeat_gtf}'\n";

fi

# check if we really need to build the STAR index

if [[ -f "${star_genome_dir}""/Genome" && -f "${star_genome_dir}""/geneInfo.tab" ]]; then

if [[ -n ${verbose} ]]; then

printf "STAR index found in '${star_genome_dir}' - using this index\n";

fi

else

# generate STAR index

# first check if genome file is .gz compressed, if so, we must uncompress it

if [[ -n ${verbose} ]]; then

printf "Beginning STAR run on '${genome_fasta_file}'\n";
printf "with '${gene_repeat_gtf}'\n";

fi

  if [[ "${genome_fasta_file}" == *".gz" ]]; then

  uncomp_genome_fasta=$(basename "${genome_fasta_file}" ".gz");
  gzip -dc ${genome_fasta_file} > ${uncomp_genome_fasta};

  "${path_to_star}"STAR --runThreadN ${starthreads} --runMode genomeGenerate --genomeDir "${star_genome_dir}" \
    --genomeFastaFiles "${uncomp_genome_fasta}" --sjdbGTFfile "${gene_repeat_gtf}";

  else

  mkdir -p "${star_genome_dir}";

  "${path_to_star}"STAR --runThreadN ${starthreads} --runMode genomeGenerate --genomeDir "${star_genome_dir}" \
    --genomeFastaFiles "${genome_fasta_file}" --sjdbGTFfile "${gene_repeat_gtf}";

  fi
fi

# create Ensembl ID to gene name list file

cat << 'ENS_SCRIPT' > ensembl_ID_to_gname.awk
# scan the gencode_gene_gtf file, extracting Ensemble IDs
# and gene names from the attributes, write a tab separated
# list of these

BEGIN{FS="\t";}

$1!~/#/&&$3=="gene"{ns = split($9,s9," ");
for (i = 1; i < ns; i+=2)
  {
  if (index(s9[i],"gene_id") > 0)
    ensid = substr(s9[i+1],2,length(s9[i+1])-3);
  if (index(s9[i],"gene_name") > 0)
    {
    printf("%s\t%s\n",ensid,substr(s9[i+1],2,length(s9[i+1])-3));
    break;
    }
  }
}
ENS_SCRIPT

if [[ -n ${verbose} ]]; then

printf "Creating Ensembl ID vs gene name file '${ensid_vs_gname}'\n";

fi

awk -f ensembl_ID_to_gname.awk "${gencode_gene_gtf}" > "${ensid_vs_gname}";

printf "genome setup complete\n";

if [[ -n delete_temp_files ]]; then

echo "Deleting scripts";

rm make_uniq_gtf_ex_allfields.awk ensembl_ID_to_gname.awk

fi

if [[ -n ${verbose} ]]; then

printf "\nRepExpress build completed\n";

fi

exit 0;

else

printf "Can't open parameter file '$1'\n";

exit 1;

fi
