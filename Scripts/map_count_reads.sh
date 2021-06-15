#!/bin/bash
#
# RepExpress script to map user reads to genome with
#  associated gtf files and run Stringtie and Featurecounts
#  on the results.
#
# requires two parameters: the name of a file containing
# basic run time information and the name of a run-specific
# file defining the mapping output directory and the read
# files.
#
# Peter Stockwell: Feb-2021
#

# check for parameter:

if [[ -z $1 || -z $2 ]]; then

printf "This script needs two parameters:\n";
printf "  1. Name of basic run parameter file\n";
printf "  2. Name of run-specific info file\n";

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
printf "  mapping and counting runs\n\n";

printf "Reading basic parameters from '$1'\n";
printf "Reading run-specific information from '$2'\n";

fi

# Do some sanity checks for necessary values

if [[ -z "${read1_fastq}" ]]; then

printf "\nError: read1_fastq file not defined in run parameter file '$2'\n";

exit 1;

fi

# check that mapping output dir exists, create if not

if [[ -z "${mapping_output_dir}" && "${mapping_output_dir}" != "." ]]; then

mkdir -p "${mapping_output_dir}";

fi

# generate a file name header for STAR output

if [[ "${read1_fastq}" == *".gz" ]]; then

r1_fq_base=$(basename "${read1_fastq}" ".fastq.gz");

else

r1_fq_base=$(basename "${read1_fastq}" ".fastq");

fi

# generate output file name

if [[ -n "${mapping_output_dir}" ]]; then

star_out_prefix="${mapping_output_dir}""/""${r1_fq_base}""_";

else

star_out_prefix="${r1_fq_base}""_";

fi

star_out_bam_name="${star_out_prefix}""Aligned.sortedByCoord.out.bam";
#star_out_bam_name="${star_out_prefix}""Aligned.out.bam";

# check if STAR has already been run, to avoid repetition

if [[ -f "${star_out_bam_name}" ]]; then

if [[ -n ${verbose} ]]; then

printf "STAR output file '${star_out_bam_name}' already exists, using this file\n"

fi

else

if [[ -n ${verbose} ]]; then

printf "running STAR, results to dir '${mapping_output_dir}'\n";

fi

# run STAR on reads:

"${path_to_star}"STAR --runThreadN "${starthreads}" --genomeDir "${star_genome_dir}" \
  --readFilesIn "${read1_fastq}" "${read2_fastq}" \
  --outFileNamePrefix "${star_out_prefix}"  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH   HI   AS   nM   XS \
  --outFilterMultimapNmax "${starfiltermax}" \
  --winAnchorMultimapNmax "${staranchormax}"  --quantMode GeneCounts;

fi

# run FeatureCounts on results, Multicount, then Unique count

if [[ -n ${verbose} ]]; then

printf "Running ${path_to_featurecounts}featureCounts on\n  '${star_out_bam_name}', multi-run\n"

fi

"${path_to_featurecounts}"featureCounts -a "${ucsc_repeats_uniq_gtf}" \
  -o "${star_out_prefix}""M_FC.txt" -t "exon" -f -p -M -O --fraction \
  -g "gene_id" ${featurecounts_overlap} ${featurecounts_strandedness} \
  -T "${featurecounts_threads}" "${star_out_bam_name}"

if [[ -n ${verbose} ]]; then

printf "Running ${path_to_featurecounts}featureCounts on\n  '${star_out_bam_name}', unique-run\n"

fi

"${path_to_featurecounts}"featureCounts -a "${ucsc_repeats_uniq_gtf}" \
  -o "${star_out_prefix}""U_FC.txt" -t "exon" -f -p -O \
  -g "gene_id" ${featurecounts_overlap} ${featurecounts_strandedness} \
  -T "${featurecounts_threads}" "${star_out_bam_name}"

# stringtie runs

if [[ -n ${verbose} ]]; then

printf "Running ${path_to_stringtie}stringtie on\n  '${star_out_bam_name}'\n"

fi

stringtie_ga_out_file="${r1_fq_base}""_gene_abund_e.out";

"${path_to_stringtie}"stringtie "${star_out_bam_name}" -G "${gencode_gene_gtf}" \
  -o "${star_out_prefix}""stringtie_e.gtf" -B -p "${stringtie_threads}" \
  -e -A "${stringtie_ga_out_file}";

# combine featureCounts Unique and Multi files, derive TPMs
# and join

cat << 'TPM_SCRIPT' > append_tpm_FC.awk
# append_tpm_FC.awk: take a FeatureCount gene abundance file
# and calculate the TPM for each line, appending the value
# to the output.  This process requires deriving a parameter
# for TPM calculations, requiring the file to be pre-scanned
# 
# usage: awk -f append_tpm_FC.awk <FeatureCount_abundance_file>
#
# command line options:
#   colhdr: header for appended column
#   nonzerotpms=1 : only write non-zero tpm values
#
# Peter Stockwell: Aug-2020

BEGIN{atotal = 0.0;
nonzerotpms = 0;
}

NR==2{printf("%s\t%s\n",$0,colhdr);
atotal = get_tpm_param(FILENAME);
if (atotal <= 0.0)
  {
  printf("TPM calculation failed, check setting of featurecounts_strandedness (-p) in mapping parameters file\n");
  exit(1);
  }
}

NR>2{if($(NF-1)+0!=0)
  {aparam=$NF/$(NF-1);
  tpm = aparam * 1.0e6 / atotal;
  }
  else
    tpm = 0.0;
if ((tpm > 0.0) || (!nonzerotpms))
  printf("%s\t%.4f\n",$0,tpm);
}

function get_tpm_param(filenam)
# scan filenam to derive tpm parameter
{
cmd = sprintf("awk -f get_tpm_parameter.awk %s",filenam);
if (cmd | getline ret > 0)
  {
  close(cmd);
  ns = split(ret,rsplit);
# printf("scanned '%s' getting param=%s\n",filenam,ret);
  return(rsplit[ns] + 0.0);
  }
else
  {
  printf("TPM parameter scan for '%s' failed\n",filenam);
  exit(1);
  }
}
TPM_SCRIPT

cat << 'SUBSCRIPT' > get_tpm_parameter.awk
# return the 'A' value for TPM from FeatureCounts output files.

BEGIN{totscaled=0.0;}

NR>2{if($(NF-1)!=0)totscaled+=$NF/$(NF-1);}

END{printf("%.2f\n",totscaled*1.0e3);}
SUBSCRIPT

# apply these scripts to featureCounts output and sort by uniq id

if [[ -n ${verbose} ]]; then

printf "Generating TPMs for Multi and Unique featureCounts runs\n"

fi

awk -f append_tpm_FC.awk colhdr="M_TPM" nonzerotpms=1 "${star_out_prefix}""M_FC.txt" | \
  sort -k 1,1 > "${star_out_prefix}""M_FC.tpm";
awk -f append_tpm_FC.awk colhdr="U_TPM" nonzerotpms=1 "${star_out_prefix}""U_FC.txt" | \
  sort -k 1,1 > "${star_out_prefix}""U_FC.tpm";

# extract a list of non_zero tpms, then join

if [[ -n ${verbose} ]]; then

printf "Combining Multi and Unique TPMs into '${star_out_prefix}M_U_FC.tpm'\n"

fi

cat "${star_out_prefix}""M_FC.tpm" "${star_out_prefix}""U_FC.tpm" | \
  cut -f 1 | sort -u | \
  join -a 1 -e "-" -o '0,2.2,2.3,2.4,2.5,2.6,2.7,2.8' - "${star_out_prefix}""M_FC.tpm" | \
  join -a 1 -e "-" -o '0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8' - "${star_out_prefix}""U_FC.tpm" | \
  tr " " "\t" > "${star_out_prefix}""M_U_FC.tpm";

# sort by tpm ('M' value) and
# put unique RE element ID at end of line,
# then run identgeneloc on this, to locate proximal genes

cat << 'SCRIPT2' > reorder_cols.awk
# put col1 to end of line

{for (i = 2; i<= NF; i++)
  printf("%s\t",$i);
printf("%s\n",$1);
}
SCRIPT2

cat << 'SCRIPT3' > add_gene_names_to_genloc.awk
# add_gene_names_to_genloc.awk: read an ENSEMBL gene_id
# from a column (16) in input and search for the corresponding
# gene name.  Append to line.
#
#

BEGIN{ensid_col = 15;
ensid_to_gene_file = "./ensid_vs_gname.txt";
while (getline ret < ensid_to_gene_file > 0)
  {
  ns = split(ret,rsplit);
  if (!(rsplit[1] in gnames))
    gnames[rsplit[1]] = rsplit[ns];
  }
close(ensid_to_gene_file);
}

$1~/#/{printf("%s\tgene\n",$0);}

$1!~/#/{ens_name = substr($ensid_col,2,length($ensid_col)-2);
if (ens_name in gnames)
  printf("%s\t%s\n",$0,gnames[ens_name]);
else
  printf("%s\t-\n",$0);
}
SCRIPT3

cat << 'APND_ST_GA' > appnd_stringtie_ga.awk
# appnd_stringie_ga.awk: script to read counts and TPMs from Stringtie gene abundance file
# into an array, indexed by the Ensemble gene ID, then append matching values to a
# file from identgeneloc with (or without) gene names appended.  This script
# will also append the gene name from Stringtie
#
# Peter Stockwell: Sep-2020
#

BEGIN{if (length(stie_file_name) <= 0)
  {
  printf("script needs a Stringtie file name as a -v stie_file_name= command option\n");
  exit(1);
  }
else
  {
  while (getline ret < stie_file_name > 0)
    {
    ns = split(ret,rsplit);
    ens_gname = rsplit[1];
    if (!(ens_gname in st_tpms))
      {
      st_tpms[ens_gname] = rsplit[ns];
      gnames[ens_gname] = rsplit[2];
      }
    }
  close(stie_file_name);
  ensid_col = 16;
  }
}

$1~/#/{printf("%s\tgene\tst_tpm\n",$0);}

$1!~/#/{ens_name = substr($ensid_col,2,length($ensid_col) - 2);
if (ens_name in st_tpms)
  printf("%s\t%s\t%s\n",$0,gnames[ens_name],st_tpms[ens_name]);
else
  printf("%s\t-\t-\n",$0);
}
APND_ST_GA

# generate a header line for output

printf "#Chromosome\tstart\tend\tstrand\tlength\thits\tM_TPM\tU_TPM\tRE_uniq_ID\tDistToGene\tLocationWRTgene\tGeneSense\tEnsemblID\tGeneCoord\tGeneName\n" > "${star_out_prefix}""M_U_FC_tpm.genloc";

# reorder columns for identgeneloc, cut to remove unwanted columns,
# then append gene name.

if [[ -n ${verbose} ]]; then

printf "Reordering columns and appending gene names,\n writing to '${star_out_prefix}M_U_FC_tpm.genloc'\n";

fi

# we need the ensid_vs_gname.txt file locally, so copy it from where it
# was made

if [[ ! -f "./ensid_vs_gname.txt" ]]; then

cp "${repeat_gene_gtf_dir}""ensid_vs_gname.txt" "./ensid_vs_gname.txt";

fi

# We won't sort this genloc file by tpm, since this makes it more complicated down
# the track for combining expression values

cat "${star_out_prefix}""M_U_FC.tpm" | awk -f reorder_cols.awk | \
"${path_to_dmap}"identgeneloc -T -f "${gencode_gene_gtf}" -i -C 9 -a "transcript" -A gene_id  -r - | \
  awk -f add_gene_names_to_genloc.awk ensid_col=15 ensid_to_gene_file="${ensid_vs_gname}" | \
  cut -f 1-11,14- >> "${star_out_prefix}""M_U_FC_tpm.genloc";

# Now sort this genloc file by col 7 (Multi featureCount tpm) and
# append stringtie gene abundance figures to produce a sorted gene abundance file

if [[ -n ${verbose} ]]; then

printf "Sorting by Multi TPM and appending stringtie gene abundance,\n"
printf " writing to '${star_out_prefix}sorted_gene_abund.genloc'\n";

fi

printf "#Chromosome\tstart\tend\tstrand\tlength\thits\tM_TPM\tU_TPM\tRE_uniq_ID\tDistToGene\tLocationWRTgene\tGeneSense\tEnsemblID\tGeneCoord\tGeneName\tStringtieAbund\n" > "${star_out_prefix}""sorted_gene_abund.genloc";

tail -n +2 "${star_out_prefix}""M_U_FC_tpm.genloc" | cut -f 1-14 | sort -nr -k 7,7  | \
  awk -v stie_file_name="${stringtie_ga_out_file}" -f appnd_stringtie_ga.awk ensid_col=13 \
  >> "${star_out_prefix}""sorted_gene_abund.genloc";

# tidy up scripts if required

if [[ -n ${delete_temp_files} ]]; then

rm append_tpm_FC.awk
rm get_tpm_parameter.awk
rm reorder_cols.awk
rm add_gene_names_to_genloc.awk
rm appnd_stringtie_ga.awk

fi

if [[ -n ${verbose} ]]; then

printf "\nMapping and counting run completed.\n\n";

printf "STAR output in '${star_out_bam_name}'\n";
printf "featureCounts multi output in '"${star_out_prefix}""_M_FC.txt"'\n";
printf "featureCounts unique output in '"${star_out_prefix}""_U_FC.txt"'\n";
printf "stringtie output in '"${star_out_prefix}""stringtie_e.gtf"'\n";
printf "featureCounts multi & unique TPMs in '"${star_out_prefix}""M_U_FC.tpm"'\n";
printf "identgeneloc output with multi TPM in '"${star_out_prefix}""M_U_FC_tpm.genloc"'\n\n";
printf "'"${star_out_prefix}""M_U_FC_tpm.genloc"' sorted with stringtie gene abundance in\n   '"${star_out_prefix}""sorted_gene_abund.genloc"'";

fi

exit 0;

else

printf "Can't open one or both parameter files '$1' or '$2'\n";

exit 1;

fi
