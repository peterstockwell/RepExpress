#!/bin/bash
#
# RepExpress script to produce a matrix of TE expression counts
# values from mapping runs performed via the map_count_reads.sh
# script, which will previously have been run on all the
# necessary samples.
#
# requires two parameters: the name of the basic run parameters
# used for the mapping and the name of a combine parameter file
# defining the information needed
#
# Peter Stockwell: Mar-2021
#

# check for parameter:

if [[ -z $1 || -z $2 ]]; then

printf "This script needs two parameters:\n";
printf "  1. Name of basic run parameter file\n";
printf "  2. Name of file parameters for combine run\n";

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
printf "sample Gene Expression hits comparison\n\n";

printf "Reading basic parameters from '$1'\n";
printf "Reading geneloc combine parameters from '$2'\n";

fi

cat << 'CT_SCRIPT' > combine_files_mk3.awk
# to combine tpm or feature files with join, filling absent
# fields.
# required input: file of sample tpm names, one
# per line.

BEGIN{nullfld = "-";
  fcount = 0;
  id_list_file = "all_TE_ids.txt";
  ostring = "0,1.2,";
  countcol = 3;
  basecmd = sprintf("join -a 1 -e \"%s\"",nullfld);
  prvcmd = sprintf("%s -o '%s2.%d' %s ",basecmd,ostring,countcol,id_list_file);
  outfile = "combined.txt";
  cmd = "";
  prvfile = "";
  hdr = "#TE_id\tSense";
}

{
ns = split($0,s0,"/");
hdr = hdr "\t" substr(s0[ns],1,length(s0[ns])-9);
ostring = "0,";
for (i = 1; i <= NR+1; i++)
  ostring = ostring sprintf("1.%d,",i+1);
cmd = cmd prvcmd;
prvcmd = sprintf(" %s | %s -o '%s2.%d' - ",$0,basecmd,ostring,countcol);
prvfile = $0;
}

END{printf("echo \"%s\" > %s\n",hdr,outfile);
printf("%s %s | tr \" \" \"\\t\" >> %s\n",cmd,prvfile,outfile);
}
CT_SCRIPT

cat << 'NAME_SCRIPT' > fname_build.awk
# take a file of sample fastq names and
# generate a file of corresponding tpm output
# files for each

BEGIN{tpm_file_location="./";}

/.fastq.gz$/{printf("%s\n",rejigname($0,"_at_r1.fastq.gz","_at_r1_M_U_FC_tpm.genloc"));}

/.fastq$/{printf("%s\n",rejigname($0,"_at_r1.fastq","_at_r1_M_U_FC_tpm.genloc"));}

/.fq$/{printf("%s\n",rejigname($0,"_at_r1.fq","_at_r1_M_U_FC_tpm.genloc"));}

/.fq.gz$/{printf("%s\n",rejigname($0,"_at_r1.fq.gz","_at_r1_M_U_FC_tpm.genloc"));}

function rejigname(oldname,oldextension,newextension)
# return the name with the new extension appended in place of old
{
ns = split(oldname,namesplit,"/");
newname=sprintf("%s/%s%s",tpm_file_location,
                 substr(namesplit[ns],1,length(namesplit[ns])-length(oldextension)),
                 newextension);
return(newname);
}
NAME_SCRIPT

cat << 'COL_SELECT' > select_genloc_cols.awk
# generate a script to select required columns from
# geneloc output files
# Note: this version uses awk

# input file is a list of files containing
# column-selected geneloc values.
# reject line 1 of input files.

BEGIN{lcnt=0;
sorted_id_file="all_TE_ids.txt";}

{ns = split($0,s0,"/");
flist[NR] = sprintf("%s_c946.txt",s0[ns]);
printf("awk 'NR>1{printf(\"%%s\\t%%s\\t%%s\\n\",$9,$4,$6);}' %s > %s\n",
        $0,flist[NR]);
lcnt = NR;
}

END{printf("cut -f 1,2");
for (i=1; i <= lcnt; i++)
  printf(" \\\n %s",flist[i]);
printf(" | sort -u -k 1,1 > %s\n",sorted_id_file);
for (i = 1; i <= lcnt; i++)
  {
  printf("rm %s\n",flist[i]) > "remove_tmp_files.sh";
  printf("%s\n",flist[i]) > "genloc_colsel_files.txt";
  }
close("remove_tmp_files.sh");
close("genloc_colsel_files.txt");
}
COL_SELECT

if [[ -n "${combine_output_dir}" ]]; then

mkdir -p "${combine_output_dir}";

fi

# now process either the fastq name list or the geneloc name list

if [[ -n ${fastq_name_file} ]]; then

  if [[ -f ${fastq_name_file} ]]; then

    if [[ -n ${verbose} ]]; then

    printf "Processing '${fastq_name_file}' of fastq files for geneloc names\n";

    fi

  if [[ -n ${verbose} ]]; then

    printf "creating list of IDs to combine - this may take a few minutes\n";

  fi

  awk -f fname_build.awk genloc_file_location="${genloc_file_dir}" "${fastq_name_file}" |  \
  awk -f select_genloc_cols.awk | /bin/sh;

  awk -f combine_files_mk3.awk outfile="${fastq_name_file}""_combined.txt" "genloc_colsel_files.txt" \
    > "${combine_output_dir}""combine_genloc.sh";

  else

  printf "Can't read fastq name file '${fastq_name_file}'\n";

  exit 1;

  fi

elif [[ -n ${genloc_name_file} ]]; then

  if [[ -f ${genloc_name_file} ]]; then

    if [[ -n ${verbose} ]]; then

    printf "Processing '${genloc_name_file}' of genloc names\n";

    fi

  if [[ -n ${verbose} ]]; then

    printf "creating list of IDs to combine - this may take a few minutes\n";

  fi

  awk -f select_genloc_cols.awk "${genloc_name_file}" | /bin/sh;

  awk -f combine_files_mk3.awk outfile="${genloc_name_file}""_combined.txt" "genloc_colsel_files.txt" \
    > "${combine_output_dir}""combine_genloc.sh";

  else

  printf "Can't read genloc name file '${genloc_name_file}'\n";
  exit 1;

  fi

fi

if [[ -n ${verbose} ]]; then

printf "Running script '"${combine_output_dir}""combine_genloc.sh"'\n\n"

fi

/bin/sh "${combine_output_dir}""combine_genloc.sh";

if [[ -n ${delete_temp_files} ]]; then

/bin/sh remove_tmp_files.sh;

rm remove_tmp_files.sh
rm combine_files_mk3.awk
rm fname_build.awk
rm select_genloc_cols.awk
rm genloc_colsel_files.txt
#rm "${combine_output_dir}""combine_genloc.sh"

fi

if [[ -n ${verbose} ]]; then

printf "\nTE expression matrix build completed.\n\n";

fi

if [[ -n ${genloc_name_file} ]]; then

printf "Results written to '"${genloc_name_file}""_combined.txt"'\n\n";

else

printf "Results written to '"${fastq_name_file}""_combined.txt"'\n\n";

fi

exit 0;

else

printf "Can't open one or both parameter files '$1' or '$2'\n";

exit 1;

fi
